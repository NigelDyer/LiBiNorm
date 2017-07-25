// ***************************************************************************
// LiBiTools.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Assorted tools associated with LiBiNorm
// ***************************************************************************

#include <fstream>
#include "api/BamReader.h"


#include "libCommon.h"
#include "stringEx.h"
#include "fastaFile.h"

#include "Options.h"
#include "FeatureFileEx.h"
#include "GeneCountData.h"
#include "LiBiTools.h"

using namespace std;
using namespace BamTools;

#define LOOKAHEAD 1000
#define MAX_READS 1000

string prefix(const string & s)
{
	return s.substr(0, s.find('.'));
}


int LiBiTools::landMain(int argc, char **argv)
{
	stringEx land_filename1, land_filename2, gff_filename, bamFileName;
	stringEx id_attribute;

	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("\nUsage: LiBiNorm land -g gff_file -b bamfile landscapefile1 landscapefile2\n\n");
		printf("This program compares two landscape files and produces a single file that combines the data from both\n");
		printf("The gff file is used to provide cooedinates for the genes so that they can be viewed easily on IGV\n");
		printf("The bam file allows the correct chromosome names to be used, the mapping being done based on chromosome length\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");
	int ni = 1;
	while (ni < argc - 2)
	{
		if (strcmp(argv[ni], "-g") == 0)
		{
			gff_filename = argv[++ni];
		}
		else if (strcmp(argv[ni], "-b") == 0)
		{
			bamFileName = argv[++ni];
		}
		else if (strcmp(argv[ni], "-i") == 0)
		{
			id_attribute = argv[++ni];
		}
		else
			exitFail("Unrecognised option", argv[ni]);
		ni++;
	}
	land_filename1 = argv[argc - 2];
	land_filename2 = argv[argc - 1];

	featureFileEx genomeDef;

	if (!id_attribute)
	{
		if (gff_filename.suffix().startsWith("gff"))
			id_attribute = DEFAULT_GFF_ID_ATTRIBUTE;
		else
			id_attribute = DEFAULT_GTF_ID_ATTRIBUTE;
	}

	stringEx feature_type = DEFAULT_FEATURE_TYPE_EXON;

	BamReader reader;
	RefVector references;

	if (!reader.Open(bamFileName))
		exitFail("Could not open input BAM files: ", bamFileName);

	// retrieve 'metadata' from BAM files.
	references = reader.GetReferenceData();
	for (auto & i : references)
	{
		if (strncasecmp(i.RefName.c_str(), "chr", 3) == 0)
			i.RefName = i.RefName.substr(3);
		genomeDef.addToChromosomeMap(i.RefLength, i.RefName);
	}

	progMessage("Reading feature file");
	if (!genomeDef.open(gff_filename, id_attribute, feature_type))
		exitFail("Could not open feature file: ", gff_filename);


	GeneCountData geneCounts1, geneCounts2,geneCountsDummy;
	progMessage("Reading first landscape file");
	geneCounts1.loadData(land_filename1);
	progMessage("Reading second landscape file");
	geneCounts2.loadData(land_filename2);

	progMessage("Indexing feature file");
	genomeDef.index(geneCountsDummy,false);

	progMessage("Outputting results");
	TsvFile resFile,missingGenesFile,extraGenesFile;
	resFile.open("combined.txt");
	missingGenesFile.open("missingGenes.txt");
	extraGenesFile.open("extraGenes.txt");

	readPositionDataClass::Iterator it1 = geneCounts1.readPositionData.begin();
//	bool used1 = geneCounts1.info[it1.readGeneAttributes().index].useForParameterEstimation;

	readPositionDataClass::Iterator it2 = geneCounts2.readPositionData.begin();
//	bool used2 = geneCounts1.info[it2.readGeneAttributes().index].useForParameterEstimation;

	while ((it1 != geneCounts1.readPositionData.end()) && (it2 != geneCounts2.readPositionData.end()))
	{
		std::map<std::string, geneData>::iterator geneIterator;

		bool match = (prefix(it1->first) == prefix(it2->first));
		if (!match)
		{
			readPositionDataClass::Iterator it1a = it1;
			for (int i = 0; (i < LOOKAHEAD) && (!match) && (++it1a != geneCounts1.readPositionData.end()); i++)
			{
				if (prefix(it1a.readGeneName()) == prefix(it2.readGeneName()))
				{
					//	There is a match for it2 later on, so it1 is a singleton.  Output it, and all entries up
					//	until the subsequent match 
					for (int j = 0; j <= i; j++)
					{
						string position("Unknown");
						string pref = prefix(it1.readGeneName());
						geneIterator = genomeDef.genes.lower_bound(pref);
						if (prefix(geneIterator->first) == pref)
							position = _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start);
						string used = geneCounts1.info[it1.readGeneAttributes().index].useForParameterEstimation ? "Y" : "";

						int size = it1.readGeneAttributes().positions[0].size() + it1.readGeneAttributes().positions[1].size();
						if (used == "Y")
						{
							missingGenesFile.print(it1.readGeneName(), _s(size, " plus"), it2.readGeneAttributes().positions[0]);
							missingGenesFile.print(it1.readGeneName(), _s(size, " minus"), it2.readGeneAttributes().positions[1]);
						}

						if (it1.readGeneAttributes().positions[0].size() > MAX_READS)
							it1.readGeneAttributes().positions[0].resize(MAX_READS);
						if (it1.readGeneAttributes().positions[1].size() > MAX_READS)
							it1.readGeneAttributes().positions[1].resize(MAX_READS);
						if (size)
						{
							resFile.print(1, position, used, it1.readGeneName(), it1.readGeneAttributes().positions[0]);
							resFile.print(1, position, used, it1.readGeneName(), it1.readGeneAttributes().positions[1]);
							resFile.print();
						}
						it1++;
					}
					match = true;
				}
			}
			if (!match)
			{
				auto it2a = it2;
				for (int i = 0; (i < LOOKAHEAD) && (!match) && (++it2a != geneCounts2.readPositionData.end()); i++)
				{
					if (prefix(it2a.readGeneName()) == prefix(it1.readGeneName()))
					{
						for (int j = 0; j <= i; j++)
						{
							string position("Unknown");
							string pref = prefix(it2.readGeneName());
							geneIterator = genomeDef.genes.lower_bound(pref);
							if (prefix(geneIterator->first) == pref)
								position = _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start);

							string used = geneCounts2.info[it2.readGeneAttributes().index].useForParameterEstimation ? "Y" : "";

							int size = it2.readGeneAttributes().positions[0].size() + it2.readGeneAttributes().positions[1].size();
							if (used == "Y")
							{
								extraGenesFile.print(it2.readGeneName(), _s(size, " plus"), it2.readGeneAttributes().positions[0]);
								extraGenesFile.print(it2.readGeneName(), _s(size, " minus"), it2.readGeneAttributes().positions[1]);
							}

							if (it2.readGeneAttributes().positions[0].size() > MAX_READS)
								it2.readGeneAttributes().positions[0].resize(MAX_READS);
							if (it2.readGeneAttributes().positions[1].size() > MAX_READS)
								it2.readGeneAttributes().positions[1].resize(MAX_READS);
							if (size)
							{
								resFile.print(2, used, position, it2.readGeneName(), it2.readGeneAttributes().positions[0]);
								resFile.print(2, used, position, it2.readGeneName(), it2.readGeneAttributes().positions[1]);
								resFile.print();
							}
							it2++;
						}
						match = true;
					}
				}
			}
		}

		string position("Unknown");
		string pref = prefix(it1.readGeneName());
		geneIterator = genomeDef.genes.lower_bound(pref);
		if (prefix(geneIterator->first) == pref)
			position = _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start);
		string used = geneCounts1.info[it1.readGeneAttributes().index].useForParameterEstimation ? "Y" : "";

		if (it1.readGeneAttributes().positions[0].size() > MAX_READS)
			it1.readGeneAttributes().positions[0].resize(MAX_READS);
		if (it1.readGeneAttributes().positions[1].size() > MAX_READS)
			it1.readGeneAttributes().positions[1].resize(MAX_READS);
		resFile.print(1, used, position, it1.readGeneName(), it1.readGeneAttributes().positions[0]);
		resFile.print(1, used, position, it1.readGeneName(), it1.readGeneAttributes().positions[1]);
		if (!match)
			resFile.print();

		position = "Unknown";
		pref = prefix(it2.readGeneName());
		geneIterator = genomeDef.genes.lower_bound(pref);
		if (prefix(geneIterator->first) == pref)
			position = _s("chr", geneIterator->second.chromosome, ":", geneIterator->second.regions[0]->start);

		used = geneCounts2.info[it2.readGeneAttributes().index].useForParameterEstimation ? "Y" : "";

		if (it2.readGeneAttributes().positions[0].size() > MAX_READS)
			it2.readGeneAttributes().positions[0].resize(MAX_READS);
		if (it2.readGeneAttributes().positions[1].size() > MAX_READS)
			it2.readGeneAttributes().positions[1].resize(MAX_READS);
		resFile.print(2, used, position, it2.readGeneName(), it2.readGeneAttributes().positions[0]);
		resFile.print(2, used, position, it2.readGeneName(), it2.readGeneAttributes().positions[1]);
		resFile.print();

		it1++;
		it2++;

	}
	return EXIT_SUCCESS;

}

int LiBiTools::landMain2(int argc, char **argv)
{
	stringEx landFilename, geneFilename;

	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("\nUsage: LiBiNorm land2 landscapefile geneFile\n\n");
		printf("This program extracts data from a landscape file\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");

	landFilename = argv[argc - 2];
	geneFilename = argv[argc - 1];

	featureFileEx genomeDef;
	stringEx id_attribute = DEFAULT_GFF_ID_ATTRIBUTE,
		feature_type = DEFAULT_FEATURE_TYPE_EXON;


	GeneCountData geneCounts;
	geneCounts.loadData(landFilename);


	ifstream file;
	file.open(geneFilename);

	if (!file.is_open())
	{
		progMessage("Unable to read gene list from ", geneFilename);
		return EXIT_FAILURE;
	}

	setEx<string> genes;

	string line, gene;
	while (!file.eof())
	{
		std::getline(file, line);
		parser(line, " \n\r", gene);
		genes.emplace(gene);
	};

	TsvFile output;
	output.open(landFilename.replaceSuffix(".subset.txt"));
	TsvFile output2;
	output2.open(landFilename.replaceSuffix(".unused.txt"));
	TsvFile genelist;
	genelist.open(landFilename.replaceSuffix(".extraGenes.txt"));

	for (size_t i = 1; i < geneCounts.info.size(); i++)
	{
#ifdef COUNT_IN_LANDSCAPE
		long count = geneCounts.counts[i];
#endif
		string & name = geneCounts.info[i].name;
		if (genes.contains(name))
		{
			long len = geneCounts.lengths[0][i];
			string c;
#ifdef COUNT_IN_LANDSCAPE
			c = _s(":", count);
#endif
			output.print(name, _s(len, c, " plus"), geneCounts.readPositionData[name].positions[0]);
			output.print(name, _s(len, c, " minus"), geneCounts.readPositionData[name].positions[1]);
		}
		else
		{
			long len = geneCounts.lengths[0][i];
			string c;
#ifdef COUNT_IN_LANDSCAPE
			c = _s(":", count);
#endif
			output2.print(name, _s(len, c, " plus"), geneCounts.readPositionData[name].positions[0]);
			output2.print(name, _s(len, c, " minus"), geneCounts.readPositionData[name].positions[1]);
			genelist.print(name);
		}

	}

	return EXIT_SUCCESS;

}


int LiBiTools::geneMain(int argc, char **argv)
{
	stringEx fastq_filename;
	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm genes fastqfile\n");
		printf("Extracts the gene names from a fasta file\n");
		return EXIT_SUCCESS;
	}
	if (argc < 2)
		exitFail("Insufficient arguments");

	fastq_filename = argv[argc - 1];

	fastaFileRead fasta;
	fasta.open(fastq_filename);

	TsvFile resFile;
	resFile.open(fastq_filename.replaceSuffix(".genes.txt"));

	bool OK = fasta.readEntry();
	while (OK)
	{
		resFile.print(fasta.NameStr());
		OK = fasta.readEntry();
	};
	return EXIT_SUCCESS;

}
