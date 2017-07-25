// ***************************************************************************
// LiBeDedup.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For removing duplicate entries in BamFiles
//	This mode is enabled using #define DEDUP_MODE in Options.h
// ***************************************************************************
#include "LiBiDedup.h"
#include <stdlib.h>
#include "api/BamReader.h"
#include "api/BamWriter.h"

#ifdef _WIN32
#include <crtdbg.h>
#endif

#include "libCommon.h"
#include "stringEx.h"
#include "Regions.h"

using namespace std;
using namespace BamTools;

int LiBiDedup::main(int argc, char **argv)
{

	minCount = 2;
	int maxGap = 3;
	bool readStart = true;

	stringEx bamFileName;

	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1],"-h") ==0) || (strcmp(argv[1],"--help")==0))))
	{
printf("Usage: LiBiNorm dedup [options] <BAM filename>\n");
printf("For each set of closely spaced reads keep just one that is aligned to the position with most reads\n");
printf("Options:\n");
printf("  -h, --help            show this help message and exit\n");
printf("  -t N, --threshold=N       minimum number of reads needed for one to be selected (3)\n");
printf("  -g N, --gap               max Gap between starts for them to be included in the same group (2)\n");
printf("  -r N, --readStart         use read start rather than match start\n");
printf("\n");
printf("Written by Nigel Dyer (nigel.dyer@warwick.ac.uk)\n");
		return EXIT_SUCCESS;
	}

	int ni = 1;

	if (argc < 1)
		exitFail("Insufficient arguments");

	while(ni < argc-1)
	{
		bool opt2 = false;

		if((strcmp(argv[ni], "-t") == 0) || (opt2 = (strncmp(argv[ni], "--threshold=",12) == 0)))
		{
			minCount = atoi(opt2?argv[ni]+12:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-g") == 0) || (opt2 = (strncmp(argv[ni], "--gap=",6) == 0)))
		{
			maxGap = atoi(opt2?argv[ni]+8:argv[++ni]);
		}
		else if((strcmp(argv[ni], "-r") == 0) || (opt2 = (strncmp(argv[ni], "--readStart",11) == 0)))
		{
			readStart = true;
		}
		else
		{
			exitFail("Invalid parameter: ",string(argv[ni]));
		}
		ni++;
	}

	bamFileName = argv[argc-1];

	BamReader reader;
	BamWriter writer;

	if ( !reader.Open(bamFileName) ) 
		exitFail("Could not open input BAM file: ",bamFileName);

	// retrieve 'metadata' from BAM files.
	RefVector references = reader.GetReferenceData();
	SamHeader header = reader.GetHeader();

	string filename(bamFileName.replaceSuffix(".dedup.bam"));
	if ( !writer.Open(filename,header,references )) 
		exitFail("Could not open output BAM file: ",bamFileName.replaceSuffix(".dedup.bam"));

#define MAXCOUNT 5000

	vector<BamAlignment> ba(MAXCOUNT+1);

	bool OK = reader.GetNextAlignment(ba[0]);

	while (OK)
	{
		int i = 0;
		size_t lastPos = ba[0].Position;
		if (readStart)
		{
			if(ba[i].CigarData[0].Type == Constants::BAM_CIGAR_SOFTCLIP_CHAR)
			{
				lastPos -=ba[i].CigarData[0].Length;
			}
		}

		size_t firstPos = lastPos;

		int thisChromosome = ba[0].RefID;
		map<size_t,vector<size_t>> instances;
		instances[0].push_back(0);

		while ((i < MAXCOUNT) && (OK = reader.GetNextAlignment(ba[++i])))
		{
			if (ba[i].RefID != thisChromosome)
				break;

			size_t pos = ba[i].Position;
			if (readStart)
			{
				if(ba[i].CigarData[0].Type == Constants::BAM_CIGAR_SOFTCLIP_CHAR)
				{
					pos -=ba[i].CigarData[0].Length;
				}
			}

			if (pos > lastPos + maxGap)
				break;
			lastPos = pos;

			if (i >= MAXCOUNT)
				i--;
			else
				instances[lastPos-firstPos].push_back(i);
		}
		if (i >= minCount)
		{
			size_t maxPos = 0;
			size_t maxPosCount = 0;
			for (auto & j:instances)
			{
				if (j.second.size() > maxPosCount)
				{
					maxPosCount = j.second.size();
					maxPos = j.first;
				}
			}

			if (instances[maxPos].size() == 1)
			{
				//only one, use it
				writer.SaveAlignment(ba[instances[maxPos][0]]);
			}
			else
			{
				size_t maxMatch = 0;
				size_t maxMatchLength = 0;
				for (auto & j:instances[maxPos])
				{

					regionList rl(ba[j]);
					size_t len = 0;
					for (auto & k:rl.data)
						len += (k.second.end- k.second.start);
					if (len > maxMatchLength)
					{
						maxMatchLength = len;
						maxMatch = j;
					}
				}
				writer.SaveAlignment(ba[maxMatch]);
			}
		}
		swap(ba[i],ba[0]);
		if (ba[0].RefID == -1)
			break;
	}

	writer.Close();
	return EXIT_SUCCESS;
}

