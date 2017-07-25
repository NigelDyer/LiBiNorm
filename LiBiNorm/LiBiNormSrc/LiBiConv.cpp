// ***************************************************************************
// LiBiConv.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For file conversion funnctions
// ***************************************************************************

#include "api/BamReader.h"
#include "libCommon.h"
#include "stringEx.h"
#include "LiBiConv.h"
#include "FeatureFileEx.h"

using namespace std;
using namespace BamTools;


int LiBiConv::main(int argc, char **argv)
{

	stringEx bamFileName, featureFileName;
	featureFileEx genomeDef;
	verbose = true;

	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("Usage: LiBiNorm conv alignment_file gff_file\n");
		printf("This program takes an alignment file in BAM format and a feature file in\n");
		printf("GFF format and creates a feature file with the chromosome names converted to the names\n");
		printf("used for the chromosomes in the bam file\n");
		return EXIT_SUCCESS;
	}
	if (argc < 3)
		exitFail("Insufficient arguments");

	bamFileName = argv[argc - 2];
	featureFileName = argv[argc - 1];

	progMessage("Converting ", featureFileName);

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

	if (!genomeDef.open(featureFileName, "","",true))
		exitFail("Could not open feature file: ", featureFileName);

	return EXIT_SUCCESS;

}

