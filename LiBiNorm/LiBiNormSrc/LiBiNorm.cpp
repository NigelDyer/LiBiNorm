// ***************************************************************************
// LiBiNorm.cpp (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// The top level code associated with "LiBiNorm model" mode
// ***************************************************************************

#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include "LiBiNorm.h"
#include "LiBiCount.h"
#include "LiBiDedup.h"
#include "LiBiConv.h"
#include "LiBiTools.h"
#include "LiBiVariation.h"
#include "MakeFastq.h"

using namespace std;

//	The top level 'main' routine.  This looks for the mode (count,model etc) and then creates an instance of the 
//	associated class and calls its 'main' routine, passing all but the initial mode parameter
int main(int argc, char **argv)
{
#ifdef _WIN32
	_DBG( _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF ));
#endif

	if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		//	No mode or help option produces a list of the modes
		printf("Usage: LiBiNorm <command> [options]\n");	
		printf("Commands:\n");	
		printf("     count            htseq-count replacement with optional bias correction\n");
		printf("     model            further bias correction analysis\n");
		printf("     conv             renames chromosomes in a .gff3 file to match those in a bam file\n");
#ifdef OUTPUT_FASTA_FILE
		printf("     refSeq           creates a fasta file with the gene reference sequences\n");
#endif
#ifdef OUTPUT_BED_DATA
		printf("     bed              creates a bed file listing exons used for analysis\n");
#endif
#ifdef INITIAL_VALUES
		printf("     variation        shows variation in Log Liklyhood with parameter\n");
#endif
#ifdef LIBITOOLS
		printf("     land             compares two different landscape files\n");
		printf("     land2            compares a landscape file and a gene list\n");
#endif
#ifdef DEDUP_MODE
		printf("     dedup            removes duplicates\n");
#endif
#ifdef MAKE_FASTQ_MODE
		printf("     makefastq        makes a fastq file from the bam file\n");
#endif
		printf("     -v, --version    version\n");
		printf("     -h, --help       this help\n");
	}
	else if (argc > 1)
	{
		string command = stringEx(argv[1]).toLower();
//		string command = (argv[1]);
		if (command == "count")
		{
			LiBiCount libiC;
			return libiC.main(argc-1,argv+1);
		}
		if (command == "model")
		{
			LiBiNorm norm;
			return norm.main(argc-1,argv+1);
		}
		if (command == "conv")
		{
			LiBiConv conv;
			return conv.main(argc - 1, argv + 1);
		}
#ifdef OUTPUT_FASTA_FILE
		if (command == "refseq")
		{
			LiBiTools tools;
			return tools.refSeqsMain(argc - 1, argv + 1);
		}
#endif
#ifdef OUTPUT_BED_DATA
		if (command == "bed")
		{
			LiBiTools tools;
			return tools.bedMain(argc - 1, argv + 1);
		}
#endif
#ifdef LIBITOOLS
		if (command == "variation")
		{
			LiBiVariation variation;
			return variation.main(argc - 1, argv + 1);
		}
		if (command == "land")
		{
			LiBiTools tools;
			return tools.landMain(argc - 1, argv + 1);
		}
		if (command == "land2")
		{
			LiBiTools tools;
			return tools.landMain2(argc - 1, argv + 1);
		}
		if (command == "genes")
		{
			LiBiTools tools;
			return tools.geneMain(argc - 1, argv + 1);
		}
#endif
#ifdef DEDUP_MODE
		if (command == "dedup")
		{
			LiBiDedup libiD;
			return libiD.main(argc - 1, argv + 1);
		}
#endif
#ifdef MAKE_FASTQ_MODE
		if (command == "makefastq")
		{
			MakeFastq makeFastq;
			return makeFastq.main(argc - 1, argv + 1);
		}
#endif
		if ((command == "--version") || (command == "-v"))
		{
			cout << "LiBiNorm version " << LIBINORM_VERSION << endl;
			return EXIT_SUCCESS;
		}
		exitFail("Invalid commmand:",command);
	}
	return EXIT_SUCCESS;
}




int LiBiNorm::main(int argc, char **argv)
{
	int Ngenes = -1;

	initClock();
	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if ((argc == 1) || ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0))))
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm model:    RNA-seq library bias normalisation using a landscape file   \n\n");
		printf("Usage: LiBiNorm model [options] landscapeFilename \n\n");
		printf("Options:\n");
		printf("  -h, --help            Show this help message and exit\n");
		helpCommon();
		printf("  -g N, --genes=N       Limit parameter discovery to using data from first N genes\n");
		printf("  -r N, --runs=N        Number of MCMC runs (", NUMBER_OF_MCMC_RUNS,")\n");
		printf("  -e N, --geneLength=N  Maximum length of transcripts used for normalisation\n");
		printf("                        parameter determination (", DEF_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION, ")\n");
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
		printf("  -k, --skip            Skip Nelder Mead parameter discovery stage. Use random initial values for MCMC\n");
		printf("  -s N, --mcmc=N        Length of each MCMC run (", NELDER_MCMC_ITERATIONS,")\n");
#else
		printf("  -s N, --mcmc=N        Length of each MCMC run (", MCMC_ITERATIONS, ")\n");
#endif
		printf("  -f, --full            Output complete set of MCMC run data\n");
		return EXIT_SUCCESS;
	}
	int ni = 1;
	while(ni < argc -1)
	{
		bool opt2 = false;
		if (commandParseCommon(ni, argv))
		{
		}
		else if ((strcmp(argv[ni], "-g") == 0) || (opt2 = (strncmp(argv[ni], "--genes=", 8) == 0)))
		{
			Ngenes = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
#ifndef USE_GROUPS_OF_GENES_FOR_DISCOVERY
			if (Ngenes < 10)
				exitFail("At least 10 genes must be specified");
#endif
		}
		else if ((strcmp(argv[ni], "-r") == 0) || (opt2 = (strncmp(argv[ni], "--runs=", 7) == 0)))
		{
			Nruns = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			if ((Nruns < 1) || (Nruns > 200))
				exitFail("-r values must lie between 1 and 200");
		}
		else if ((strcmp(argv[ni], "-s") == 0) || (opt2 = (strncmp(argv[ni], "--mcmc=", 7) == 0)))
		{
			Nsimu = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			if ((Nsimu < 1) || (Nsimu > 10000))
				exitFail("-s values must lie between 1 and 10000");
		}
		else if ((strcmp(argv[ni], "-e") == 0) || (opt2 = (strncmp(argv[ni], "--geneLength=", 13) == 0)))
		{
			maxGeneLength = atoi(opt2 ? argv[ni] + 13 : argv[++ni]);
		}
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
		else if ((strcmp(argv[ni], "-k") == 0) || (opt2 = (strncmp(argv[ni], "--skip", 6) == 0)))
		{
			nelderMead = false;
		}
#endif
		else
		{
			exitFail("Invalid parameter: ",argv[ni]);
		}
		ni++;
	}

	landscapeFilename = argv[argc - 1];

	//	If we specifiy the model then run the other models just once 
	Nmodels = allModels().size();
	if (theModel == none)
	{
		normalise = false;
		Nmodels = 0;
	}
#ifdef SELECT_BY_LL
	else if (theModel == noModelSpecified)
		NrunsOtherModels = Nruns;
	else if (theModel == findBestModel)
		NrunsOtherModels = Nruns;
	else
	{
//		NrunsOtherModels = (Nruns == 1) ? 0 : 1;
		NrunsOtherModels = 0;
		Nmodels = 1;
	}
#else
	else if (calcAllModels)
	{
		NrunsOtherModels = Nruns;
	}


#endif
#ifdef USE_GROUPS_OF_GENES_FOR_DISCOVERY
	geneCounts.loadData(landscapeFilename, 1000000);

	bool useGene = true;
	int geneCount = 0;
	if (Ngenes < 0)
	{
		useGene = false;
		Ngenes = -Ngenes;
	}
	for (int i = 0; i < geneCounts.info.size(); i++)
	{
		if (!useGene)
		{
			geneCounts.info[i].name = "XXX";
			geneCounts.info[i].useForParameterEstimation = false;
		}
		if (++geneCount == Ngenes)
		{
			useGene = !useGene;
			geneCount = 0;
		}
	}
#else
	allGeneCounts.loadData(landscapeFilename, Ngenes);
#endif



	if (outputFileroot)
		allGeneCounts.outputHeatmapData(outputFileroot.replaceSuffix("_bias.txt"));

	//	Load up the initial values
	if (parameterFilename)
	{
		SetInitialParamsFromFile(parameterFilename);
	}

	if (normalise)
	{

		coreParameterEstimation();

		//	If we have explicitly specified the model then use it instead
#ifdef SELECT_BY_LL
		if ((theModel == noModelSpecified) || (theModel == findBestModel))
		{
			bestModel = getBestModel();
			progMessage("Best model is ", bestModel);
			theModel = bestModel;
		}
		else
#endif
		{
			progMessage("Model selected by command line is ", theModel);
		}

		getBias(theModel, bestResults[theModel].params[logValue], allGeneCounts.lengths[0], allGeneCounts.bias);

		//	If a file root was specified then output the detailed output files
		if (outputFileroot)
		{
			string filename = outputFileroot.replaceSuffix("_expression.txt");
			if (!allGeneCounts.outputGeneCounts(filename, outputFPKM ? 3 : 2, conv(theModel)))
				exitFail("Unable to output counts to :", filename);

			printResults();
			printBias();
#ifdef PRINT_DISTRIBUTION
			printDistribution(allGeneCounts);
#endif
		}

		//	Optional print out of all of the data for the full set of mcmc runs for each model
#ifdef PRINT_MCMC_RUN_DATA
		printAllMcmcRunData();
		printConsolidatedMcmcRunData();
#endif
	}
#ifdef PRINT_DISTRIBUTION
	else if(outputFileroot)
		printDistribution(allGeneCounts);
#endif

	//	The basic count data in htseq-count format
	if (!allGeneCounts.outputGeneCounts(countsFilename, normalise?1:0, conv(theModel)))
		exitFail("Unable to output counts to :", countsFilename);

	progMessage("Data modelled");
	elapsedTime();

	if (pauseAtEnd)
	{
		string test;
		cin >> test;
	}

	return EXIT_SUCCESS;
}

