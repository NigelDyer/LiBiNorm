// ***************************************************************************
// MakeFastq.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For making Fastq Files with various optional artefacts from BamFiles
// ***************************************************************************

#include <random>
#include <chrono>
#include "MakeFastq.h"
#include "fastaFile.h"

#ifdef _DEBUG
#define BAMCACHESIZE 100
#else
#define BAMCACHESIZE 5000
#endif

unsigned seed()
{
		auto now_us = chrono::time_point_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now());
		return chrono::duration_cast<chrono::microseconds>(now_us.time_since_epoch()).count();
};


default_random_engine generator(seed());
uniform_real_distribution<double> distribution(0.0,1.0);

bamRead::bamRead(const BamAlignment & ba) 
{
	operator = (ba);
}

bamRead & bamRead::operator = (const BamAlignment & ba)
{
	name = ba.Name;
	if (ba.IsMapped() && ba.IsReverseStrand())
	{
		readSeq = sequence(ba.QueryBases).complement();
		qualData = sequence(ba.Qualities).reversed();
	}
	else
	{
		readSeq = ba.QueryBases;
		qualData = ba.Qualities;
	}
	return This;
}

void bamRead::setName(const BamAlignment & ba)
{
	name = ba.Name;
}

void bamRead::setValues(const bamRead & br)
{
	readSeq = br.readSeq;
	qualData= br.qualData;
}

//	Add a SNP by setting a random position in the read to the value at another random position
void bamRead::addSNP(double errorRate)
{
	double randVal = distribution(generator);
	if (randVal > (errorRate*readSeq.size()))
		return;

	size_t p1 = ceil(distribution(generator) * readSeq.size()-1);
	size_t p2 = ceil(distribution(generator) * readSeq.size()-1);

	readSeq[p1] = readSeq[p2];
}

//	Output a line of fastq data
void bamRead::output(FILE * f)
{
	fprintf(f,"@%s\n%s\n+\n%s\n",name.c_str(),readSeq.c_str(),qualData.c_str());
}

void bamReadCache::clear()
{
	for (auto & i: This)
		i.resetNames();
}

//	Get next alignment.  If we reach the end then go back to the beginning
inline bool GNA(BamAlignment & ba,BamReader & br)
{
	//	Returns false if we read the end and have to rewind (oe some other fault occurs), otherwise true;
	int NH;
	bool OK;
	do
	{
		OK = br.GetNextAlignment(ba);
		if (!OK)
		{
			cerr << "Rewinding" << endl;
			br.Rewind();
			br.GetNextAlignment(ba);
		}
		// While this can be rearranged, this expression states what is required most clearly 
	} while (!(!ba.GetTag("NH",NH)  || NH <= 1 || ba.IsPrimaryAlignment()));

	return OK;
}

bool MakeFastq::getNextAlignment()
{
	return GNA(ba,reader);
}

bool MakeFastq::getNextForeignAlignment()
{
	return GNA(foreignBa,foreignReader);
}

bool MakeFastq::getNextAlignmentCore()
{
	//	Cannot check for tags as these are not decoded in GetNextAlignmentCore
	return reader.GetNextAlignmentCore(ba);
}

void MakeFastq::incrementCount()
{
	if ((++count % 10000) == 0)
		cerr << "Record " << count << endl;
}

//
//	Not really part of LiBiNorm but ended up being incorporated as an option
//	This can be used to generate fastq files from BAM files with additional artefacts
int MakeFastq::main(int argc, char **argv)
{

	int size=1500000;
		
	double	mitoRate=1,unmappedRate=1,mappedRate=1,
		foreignBamDataRate = 0;;

	int	overamplified = -1;
	double overAmplifiedError = 0.002;

	size_t start = 0;

	stringEx bamFileName,outputFileRoot,foreignBamData;

	stringEx skipChromosome;

	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1],"-h") ==0) || (strcmp(argv[1],"--help")==0))))
	{
		printf("Usage: LiBiNorm makefastq [options] alignment_file\n");
		printf("Options:\n");
		printf("  -h, --help				show this help message and exit\n");
		printf("  -o <nam>, --output=<name>	The root name for the generated fastq files\n");
		printf("  -s N, --size=N			Number of reads in each fastq file in  (1500000)\n");
		printf("  -k ab, --skip=ab			skip chromosomes begining in ab\n");
		printf("  -i f, --mito=f			Proportion of mitochondrial reads from original\n");
		printf("							that are included in output (1.0)\n");
		printf("  -u f, --unmapped=f		Proportion of completely unmapped reads in output (1.0)\n");
		printf("  -m f, --mapped=f			Proportion of mapped reads in output (1.0)\n");
		printf("  -v N f, --overamplified=N,f	degree of overamplification (N) and error rate in copies (g,0.002) \n");
		printf("  -f <filename> f g, --foreign<filename>,f,g\n");
		printf("							f = Proportion of foreign DNA from bamfile <filename> included in output (0.0)\n");
		printf("							g = Error rate of duplicated reads (0.002)\n");
		printf("\n");
		printf("Written by Nigel Dyer (nigel.dyer@warwick.ac.uk)\n");
		return EXIT_SUCCESS;
	}

	FILE * f1out,* f2out;

	int ni = 1;
	bool pendingOveramplifed = false;

	if (argc < 1)
		exitFail("Insufficient arguments");

	while(ni < argc-1)
	{
		bool opt2 = false;

		if((strcmp(argv[ni], "-l") == 0) || (opt2 = (strncmp(argv[ni], "--length=",9) == 0)))
		{
			size = atoi(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-o") == 0) || (opt2 = (strncmp(argv[ni], "--output=",9) == 0)))
		{
			outputFileRoot = opt2?argv[ni]+9:argv[++ni];
		}
		else if((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--start=",8) == 0)))
		{
			start = atoi(opt2?argv[ni]+8:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-i") == 0) || (opt2 = (strncmp(argv[ni], "--mito=",7) == 0)))
		{
			mitoRate = atof(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-u") == 0) || (opt2 = (strncmp(argv[ni], "--unmapped=",11) == 0)))
		{
			unmappedRate = atof(opt2?argv[ni]+11:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-m") == 0) || (opt2 = (strncmp(argv[ni], "--mapped=",9) == 0)))
		{
			mappedRate = atof(opt2?argv[ni]+9:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-v") == 0) || (opt2 = (strncmp(argv[ni], "--overamplified=",9) == 0)))
		{
			overamplified = atoi(opt2?argv[ni]+9:argv[++ni]);
			if (opt2)
			{
				const char * comma = strchr(argv[ni]+9,',');
				if (comma)
					overAmplifiedError = atof(comma+1);
			}
			else 
				pendingOveramplifed = true;
		}
		else if((strcmp(argv[ni], "-k") == 0) || (opt2 = (strncmp(argv[ni], "--skip=",6) == 0)))
		{
			skipChromosome = opt2?argv[ni]+7:argv[++ni];
		}
		else if((strcmp(argv[ni], "-f") == 0) || (opt2 = (strncmp(argv[ni], "--foreign=",10) == 0)))
		{
			foreignBamData = opt2?argv[ni]+10:argv[++ni];
			foreignBamDataRate = atof(opt2?argv[ni]+11+foreignBamData.size():argv[++ni]);
		}
		else if (pendingOveramplifed)
		{
			overAmplifiedError = atof(argv[ni]);
		}
		else
		{
			exitFail("Invalid parameter: ",string(argv[ni]));
		}
		ni++;
	}

	bamFileName = argv[argc-1];

	if (!outputFileRoot)
		outputFileRoot = bamFileName.removeSuffix();

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM file: ",bamFileName);
	// retrieve 'metadata' from BAM files.
	RefVector references = reader.GetReferenceData();

	if (foreignBamData)
	{
		if ( !foreignReader.Open(foreignBamData) ) 
			exitFail("Could not open input BAM file: ",foreignBamData);
		for (size_t i = 0;i < start;i++)
			getNextForeignAlignment();
	}

	int mitoRef = -1;
	for (size_t i = 0;i < references.size();i++)
	{
		if ((references[i].RefName == "MT") || (references[i].RefName == "chrMT"))
		{
			mitoRef = i;
			break;
		}
	}
	if (mitoRef == -1)
		exitFail("No mitochondrial gene found");


	cerr << "Skipping reads" << endl;
	for (size_t i = 0;i < start;i++)
		getNextAlignmentCore();
	bool OK = getNextAlignment();


	//	Use binary output format so that we get unix style/n line feeds
	f1out = fopen((outputFileRoot + "_R1_001.fastq").c_str(),"wb");
	f2out = fopen((outputFileRoot + "_R2_001.fastq").c_str(),"wb");


	cerr << "Creating fastq files" << endl;
	bamReadCache oaBuffer(BAMCACHESIZE);
	size_t overAmplifyRounds = 0;
	pairedBamReads  readPair,readForeignPair;
	size_t buffIndex = 0;
	count = 0;
	while (OK && (count < size))
	{

		bool useThis = true;

		if (skipChromosome && (ba.RefID != -1) && (stringEx(references[ba.RefID].RefName).startsWith(skipChromosome)))
		{
			useThis = false;
		}

		if (ba.IsFirstMate())
		{
			double r2 = distribution(generator);

			if (!ba.IsMapped() && !ba.IsMateMapped())
			{
				if (r2 >= unmappedRate)
					useThis = false;
			}
			else if ((ba.IsMapped() && (ba.RefID == mitoRef)) || (ba.IsMateMapped() && (ba.MateRefID == mitoRef)))
			{
				if (r2 >= mitoRate)
					useThis = false;
			}
			else
			{
				if (r2 >= mappedRate)
					useThis = false;
			}
		}

		if (useThis)
		{
			if (overamplified > 1)
			{
				ba.BuildCharData();
				oaBuffer[buffIndex].addRead(ba);

				if (oaBuffer[buffIndex].namesMatch())
						buffIndex++;

				if (buffIndex >= BAMCACHESIZE)
				{
					size_t i = 0;

					size_t amp = 1;
					switch (overAmplifyRounds++ % 6)
					{
					case 0: amp = 2; break;
					default:amp = overamplified; break;
					}


					while ((i++ < (BAMCACHESIZE * amp)) && (count < size) && OK)
					{
						incrementCount();

//						int NH = -1;
						// Find a newname for the overamplified record.  Need a logic so that we only pick one 
						// read with any given name	
						do {
							OK = getNextAlignment();
						}
						while (!ba.IsFirstMate() && OK);

						//	Pick one of the records, rename it with the name from another record, add some errors and output it
						size_t pos = distribution(generator) * (BAMCACHESIZE-1);

						pairedBamReads pbr(oaBuffer[pos]);
						pbr.setName(ba);
						pbr.addSNP(overAmplifiedError);
						pbr.output(f1out,f2out);

					}
					oaBuffer.clear();
					buffIndex = 0;
				}

			}
			else
			{
				readPair.addRead(ba);

				if (readPair.namesMatch())
				{
					incrementCount();

					//	A fraction r2 of the records are replaced with reads from a selected foreign bamfile
					double r2 = distribution(generator);
					if ((foreignBamData) && (r2 < foreignBamDataRate))
					{
						while (!readForeignPair.namesMatch())
						{
							readForeignPair.addRead(foreignBa);
							getNextForeignAlignment();
						}
						readPair.setValues(readForeignPair);
						readForeignPair.resetNames();
					}
					readPair.output(f1out,f2out);
					readPair.resetNames();

				}
			}
		}

		OK = getNextAlignment();
	}

	fclose(f1out);
	fclose(f2out);

	return EXIT_SUCCESS;
}
