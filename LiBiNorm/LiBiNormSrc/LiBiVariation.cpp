// ***************************************************************************
// LiBiVariation.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// The top level code associated with "LiBiNorm count" modes
// ***************************************************************************

#include "parser.h"
#include "LiBiVariation.h"

using namespace std;

//
//	This function is only available if LIBITOOLS is defined in Options.h
int LiBiVariation::main(int argc, char **argv)
{
	stringEx landscapeFilename;
	initClock();
	if (argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm variation:  Test of effect of parameter variation\n\n");
		printf("Usage: LiBiNorm variation [Options] -N <resultsFileroot> -p <parameterFile> <landscapeFile>\n\n");
		printf("Options:\n");
		printf("  -e N, --seed=N       Ensures a specific set of reads are selected\n");
		printf("  -h, --help            show this help message and exit\n");
		helpCommon();
		return EXIT_SUCCESS;
	}
	int ni = 1;
	while (ni < argc-1)
	{
		bool opt2 = false;
		if (commandParseCommon(ni, argc, argv))
		{
		}
		else if ((strcmp(argv[ni], "-e") == 0) || (opt2 = (strncmp(argv[ni], "--seed=", 7) == 0)))
		{
			int seed = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			intRandClass::instance().reseed(seed);
		}
		else
		{
			exitFail("Invalid parameter: ", argv[ni]);
		}
		ni++;
	}

	landscapeFilename = argv[argc - 1];

	if (!parameterFilename)
		exitFail("Parameter filename must be defined");

	if (!outputFileroot)
		exitFail("-N parameter must be specified");

	//	Load up the initial values
	SetInitialParamsFromFile(parameterFilename);

	//	Load up the landscape file
	geneCounts.loadData(landscapeFilename, -1);
	geneCounts.remove_invalid_values();
	geneCounts.transferTo(geneData, MAX_READS_GENE, maxReads, maxGeneLength);

	string filename(outputFileroot.replaceSuffix("_variation.txt"));
	TsvFile mcmcResult;
	//	Now output a table with the end points of each of the chains.
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	First headers up to and including the maximum model that is run.   Always leave space
	//	for the intermediate models so the layout of the results is consistent
	mcmcResult.printStart("");
	for (modelType m : allModels())
	{
		mcmcResult.printMiddle(m);
		mcmcResult.printRepeat(headers[m].size() + 1);
	}
	mcmcResult.printEnd();

	mcmcResult.printStart("");
	for (modelType m : allModels())
		mcmcResult.printMiddle(headers[m], "chain", "");
	mcmcResult.printEnd();


	mcmcResult.printStart("Log Opt");
	for (modelType m : allModels())
		mcmcResult.printMiddle(initialValues[m], "" , "");
	mcmcResult.printEnd();

	optionsType options;

	map<modelType, paramDescriptionSet> params;
	for (modelType m : allModels())
		params[m] = GetModelParams(m);

	for (size_t p = 0; p < 2; p++)
	{
		for (double i = 0; i < 3; i += 0.1)
		{
			mcmcResult.printStart("");
			for (modelType m : allModels())
			{
				dataVec values = initialValues[m];
				values[p] = params[m][p].min + i;
				setSSfun(options, m);
				double ss1 = options.ssfun(values, geneData);
				mcmcResult.printMiddle(values, ss1, "");
			}
			mcmcResult.printEnd();
		}
		mcmcResult.print("");
	}
	return EXIT_SUCCESS;
}

