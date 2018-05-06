// ***************************************************************************
// LiBiNormCore.cpp (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// Common code associated with all of the LiBiNorm modes
// ***************************************************************************
#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <mutex>
#include <thread>
#include "LiBiNormCore.h"
#include "LiBiOptimiser.h"

using namespace std;

void LiBiNormCore::helpCommon()
{
	printf("  -n M, --normModel=M   Specifies that model M should be used rather than the default\n");
	printf("                        Model BD. M options: A,B,C,D or polyA,\n");
	printf("                        E or random,BD or smart");
#ifdef SELECT_BY_LL
	printf(",best.  best causes all models to\n");
	printf("                        be evaluated and the best, based on liklihood selected\n");
#endif
	printf("\n");
	printf("  -u FILEROOT, --normFileroot=FILEROOT\n");
	printf("                        All output summary info is sent to files with root FILEROOT\n");
	printf("  -p N, --threads=N     Number of threads for normalisation parameter\n");
	printf("                        determination (", DEF_THREADS, ")\n");
	printf("  -d N, --reads=N       Maximum number of reads using for normalisation\n");
	printf("                        parameter determination (", DEF_MAX_READS_FOR_PARAM_ESTIMATION, ")\n");
#ifndef SELECT_BY_LL
	printf("  -f, --full            Calculate full set of models\n");
#endif
#ifdef INITIAL_VALUES
	printf("  -v <filename>, --initial=FILENAME\n");
	printf("                        Set initial values for parameter discoverey from file\n");
#endif
	printf("  -q, --quiet           suppress progress report\n");
#ifdef OUTPUT_DEBUG_MESSAGES
	printf("  -w, --debug           output debug messages\n");
#endif
#ifdef PAUSE_AT_END_OPTION
	printf("  -x	                pause at end rather than simply exiting\n");
#endif
#ifdef SELECT_READ_SEED
	printf("  -y N, --seed=N        Set seed used for selecting subset of reads\n");
#endif
	printf("  -c FILENAME, --counts=FILENAME\n");
	printf("                        Name of output file. default: writes to stdout\n");

}

void LiBiNormCore::featureAndIdAttributeHelp()
{
	printf("  -t FEATURETYPE, --type=FEATURETYPE\n");
	printf("                        feature type (3rd column in GFF file) to be used, all\n");
	printf("                        features of other type are ignored (default for \n");
	printf("                        Ensemble GTF and GFF files: ", DEFAULT_FEATURE_TYPE_EXON, ")\n");
	printf("  -i IDATTR, --idattr=IDATTR\n");
	printf("                        GFF attribute to be used to identify the feature ID\n");
	printf("                        (default for Ensembl GTF files: ", DEFAULT_GTF_ID_ATTRIBUTE, "\n");
	printf("                        default for GFF3 files: ", DEFAULT_GFF_ID_ATTRIBUTE, ")\n");
}


bool LiBiNormCore::commandParseCommon(int & ni, char **argv)
{
		bool opt2 = false;
		if ((strcmp(argv[ni], "-u") == 0) || (opt2 = (strncmp(argv[ni], "--outputFileroot=", 17) == 0)))
		{
			outputFileroot = opt2 ? argv[++ni] + 17 : argv[++ni];
			return true;
		}
		if ((strcmp(argv[ni], "-n") == 0) || (opt2 = (strncmp(argv[ni], "--normModel=", 12) == 0)))
		{
			theModel = modelFromString(opt2 ? argv[++ni] + 12 : argv[++ni]);
			return true;
		}
		if ((strcmp(argv[ni], "-p") == 0) || (opt2 = (strncmp(argv[ni], "--threads=", 10) == 0)))
		{
			Nthreads = atoi(opt2 ? argv[ni] + 10 : argv[++ni]);
			if (Nthreads < 1)
				exitFail("At least 1 thread must be specified");
			return true;
		}
		if ((strcmp(argv[ni], "-d") == 0) || (opt2 = (strncmp(argv[ni], "--reads=", 8) == 0)))
		{
			maxReads = atoi(opt2 ? argv[ni] + 8 : argv[++ni]);
			if (maxReads < 1000)
				exitFail("At least 1000 reads must be specified");
			return true;
		}
		if ((strcmp(argv[ni], "-c") == 0) || (opt2 = (strncmp(argv[ni], "--counts=", 9) == 0)))
		{
			countsFilename = opt2 ? argv[ni] + 9 : argv[++ni];
			return true;
		}
		if ((strcmp(argv[ni], "-j") == 0) || (opt2 = (strncmp(argv[ni], "--fkpm", 6) == 0)))
		{
			outputFPKM = true;
			return true;
		}

#ifndef SELECT_BY_LL
		if ((strcmp(argv[ni], "-f") == 0) || (opt2 = (strncmp(argv[ni], "--full", 6) == 0)))
		{
			calcAllModels = true;
			return true;
		}
#endif
#ifdef OUTPUT_DEBUG_MESSAGES
		if ((strcmp(argv[ni], "-w") == 0) || (opt2 = (strncmp(argv[ni], "--debug", 7) == 0)))
		{
			debugPrint = true;
			return true;
		}
#endif
#ifdef PAUSE_AT_END_OPTION
		if (strcmp(argv[ni], "-x") == 0)
		{
			pauseAtEnd = true;
			return true;
		}
#endif

		if ((strcmp(argv[ni], "-q") == 0) || (opt2 = (strncmp(argv[ni], "--quiet", 7) == 0)))
		{
			verbose = false;
			return true;
		}
#ifdef INITIAL_VALUES
		if ((strcmp(argv[ni], "-v") == 0) || (opt2 = (strncmp(argv[ni], "--intial=", 9) == 0)))
		{
			parameterFilename = (opt2 ? argv[ni] + 8 : argv[++ni]);
			return true;
		}
#endif
#ifdef SELECT_READ_SEED
		if ((strcmp(argv[ni], "-y") == 0) || (opt2 = (strncmp(argv[ni], "--seed=", 7) == 0)))
		{
			int seed = atoi(opt2 ? argv[ni] + 7 : argv[++ni]);
			intRandClass::instance().reseed(seed);
			return true;
		}
#endif

		return false;
}


bool LiBiNormCore::commandParseIdAndType(int & ni, char **argv)
{
	bool opt2 = false;
	if ((strcmp(argv[ni], "-t") == 0) || (opt2 = (strncmp(argv[ni], "--type=", 7) == 0)))
	{
		feature_type = opt2 ? argv[ni] + 7 : argv[++ni];
		return true;
	}
	if ((strcmp(argv[ni], "-i") == 0) || (opt2 = (strncmp(argv[ni], "--idattr=", 9) == 0)))
	{
		id_attribute = opt2 ? argv[ni] + 9 : argv[++ni];
		return true;
	}
	return false;
}

void LiBiNormCore::SetInitialParamsFromFile(const string & filename)
{
	parseTsvFile paramFile;
	if (!paramFile.open(filename))
		exitFail("Unable to read parameters from ", filename);
	paramFile.read(initialValues);


	cout << "Need to check that SetIntialValuesFromFile works correctly" << endl;

	//	Get rid of spurious values (possibly as a result of trailing tabs in the text file)
	for (modelType m : allModels())
		initialValues[m].resize(headers[m].size());
}

//	Multiple instances of this are called, each works through the requested mcmc runs as 
//	listed in threadLoopCounts which has a pair of integers associated with each 
void LiBiNormCore::mcmcThread(optionsType options)
{
	static mutex mtx1;

	static vector<mutex> initValuesMutexes(allModels().size() + 1);
	if (nelderMead)
	{
		modelType m;
		bool modelsLeftToDo = true;
		while (modelsLeftToDo)
		{
			{
				lock_guard<mutex> lock(mtx1);

				while ((modelsLeftToDo = (nelderMeadCounter < allModels().size())))
				{
//					m = allModels()[allModels().size() - nelderMeadCounter++ -1];
					m = allModels()[nelderMeadCounter++];
					if (threadLoopCounts[m].requested > 0)
						break;
				}
				if (!modelsLeftToDo)
					break;
				initValuesMutexes[m].lock();
				progMessage("Starting intial values ", m);
			}

			LiBiOptimiser  optimiser(allGeneData, m, bestResults[m].nelderMeadIterations);
			setSSfun(options, m);
			initialValues[m] = optimiser.getParams(m, options, initialValues[m]);

			initValuesMutexes[m].unlock();

			lock_guard<mutex> lock(mtx1);
			progMessage("Finishing intial values ", m);
		}
	}

	map<modelType, loop_counts >::iterator model_iterator = threadLoopCounts.begin();
	modelType currentModel;
	mcmcRunId mcmcRun;
	while (true)
	{
		//  Look for an iteration of a loop that is yet to be done.
		{
			lock_guard<mutex> lock(mtx1);
			//	Check to see if we have done all the runs associated with this model
			//	if so then move on to the next
			while (model_iterator->second.counter > model_iterator->second.requested)
			{
				if (++model_iterator == threadLoopCounts.end())
					//	All models have been done
					return;
			}
			mcmcRun = model_iterator->second.counter++;

			currentModel = model_iterator->first;
		}
		if (nelderMead)
		{
			//	Stop here if necessary to wait for the initial values to be completed.  
			//	Only likely to happen if there are more than 6 threads
			lock_guard<mutex> lock(initValuesMutexes[currentModel]);
		}
		{
			lock_guard<mutex> lock(mtx1);
			progMessage("Starting ", currentModel, ", iteration:", mcmcRun);
		}
		setSSfun(options, currentModel);

		paramDescriptionSet params = GetModelParams(currentModel, &initialValues[currentModel], options.jumpSize);
		options.qcov = dataVec(params.size(), options.jumpSize);

		mcmc mcmcEngine;
		mcmcEngine.mcmcrun(allGeneData, params, options);

		//  Make sure only one thread at a time is outputting results
		lock_guard<mutex> lock(mtx1);

		progMessage("Finishing ", currentModel, ", iteration:", mcmcRun);

#ifdef STORE_ENDPOINTS
		Chain[currentModel].emplace(loop, mcmcEngine.chain().back());
		SSChain[currentModel].emplace(loop, mcmcEngine.sschain().back());
#endif

		//	Always store full set of results as these are needed to calculate the optimal parameters
		//	emplace/move them for efficiency
		fullResultChain[currentModel].emplace(mcmcRun, move(mcmcEngine._chain));
		fullResultSSChain[currentModel].emplace(mcmcRun, move(mcmcEngine._sschain));
	}
}

void LiBiNormCore::checkFeatureAndIdAttribute()
{
	if (!feature_type)
	feature_type = DEFAULT_FEATURE_TYPE_EXON;

	if (!id_attribute)
	{
		if (featureFileName.suffix() == "gtf")
			id_attribute = DEFAULT_GTF_ID_ATTRIBUTE;
		else if (featureFileName.suffix().startsWith("gff"))
			id_attribute = DEFAULT_GFF_ID_ATTRIBUTE;
		else
			exitFail("Unable to identify feature file type in order to specifiy default id attribute");
	}
}

bool LiBiNormCore::coreParameterEstimation()
{
	allGeneCounts.remove_invalid_values();
	allGeneCounts.transferTo(allGeneData, MAX_READS_GENE, maxReads, maxGeneLength);

	elapsedTime("Data loaded");

	optionsType options;

	options.jumpSize = MCMC_JUMP_SIZE;
	options.nsimu = Nsimu;
	options.Nruns = Nruns;

	options.sigma2 = MCMC_SIGMA;

	//	Set the number of iterations required of each of the models.
	for (modelType m : allModels())
	{
		threadLoopCounts[m].counter = 1;
		if (m == theModel)
			threadLoopCounts[m].requested = options.Nruns;
		else
			threadLoopCounts[m].requested = NrunsOtherModels;
	}

#ifdef FIXED_RESULTS
	theModel = M_FIXED_RESULTS;
	bestResults[theModel].params = dataVec{ FIXED_RESULTS };
#else

	nelderMeadCounter = 0;
	//	And then set the threads running
	if (Nthreads == 1)
	{
		//	Dont make additional threads if single threaded.  Makes debugging easier
		mcmcThread(options);
	}
	else
	{
		vector<thread> threads;
		for (size_t i = 0; i < Nthreads; i++)
			threads.emplace_back(&LiBiNormCore::mcmcThread, this, options);

		for (auto & i : threads)
			i.join();
	}

#ifdef XXXX
	Unfinished code for
		vector<int> geneL{ 500,1000,2000,4000,8000 };
	dataVec TotalReads;

	for (size_t i = 0; i < geneL.size(); i++)
	{
		for (size_t j = 0; j < transData.size(); j++)
		{
			if (abs(transData[j].length - geneL[i]) < (0.1 * geneL[i]))
				TotalReads.append(transData[j].counts[0] / transData[j].length);
		}
	}
#endif

	//********************************************************************************************
	//	Find the optimal parameter values, which are associated with results found in the last
	//	iterations of all of the MCMC runs.
	for (modelType m : allModels())
	{
		multimap <double, dataVec *> & orderedResults = allOrderedResults[m];
		if (fullResultSSChain[m].size())
		{
			bestResult & br = bestResults[m];
			size_t Nparams = headers[m].size();
			for (size_t i = 0; i < 2; i++)
				br.params[i].resize(Nparams);
			for (size_t i = 0; i < 4; i++)
				br.param_dev[i].resize(Nparams);

			//	Use the second half of teh chain to find parameters and their variation
			for (mcmcRunId i = 1; i <= fullResultSSChain[m].size(); i++)
			{
				for (size_t j = fullResultSSChain[m][i].size() / 2; j < fullResultSSChain[m][i].size(); j++)
				{
					orderedResults.emplace(fullResultSSChain[m][i][j], &fullResultChain[m][i][j]);
				}
			}

#ifdef USE_PARAMS_FROM_LOWEST_LL
			//	Find the best LL.
			br.LLresult = orderedResults.begin()->first;
			//	Get the parameters from the most likly sample
			for (size_t p = 0; p < Nparams; p++)
			{
				VEC_DATA_TYPE v = (*orderedResults.begin()->second)[p];
				br.params[logValue][p] = v;
				if (p < 4)
					br.params[absValue][p] = pow(10, v);
				else
					br.params[absValue][p] = v;
			}
#else
			{
				//	Put the Log liklyhoods in order and find the median
				auto i = orderedResults.begin();
				size_t n = 0;
				for (; n < orderedResults.size() / 2; i++, n++) {};
				br.LLresult = i->first;
			}
			{
				//	Now put each of the params in order
				auto j = orderedResults.begin();
				vector<orderedVec> orderedParams(Nparams);
				for (; j != orderedResults.end(); j++)
				{
					for (size_t p = 0; p < Nparams; p++)
						orderedParams[p].add((*j->second)[p]);
				}
				//	and then find the median
				for (size_t p = 0; p < Nparams; p++)
				{
					VEC_DATA_TYPE v = orderedParams[p].median();
					br.params[logValue][p] = v;
					if (p == 1)
#ifdef ABS_H_PARAM
						br.params[absValue][1] = v;
#else
						br.params[absValue][1] = pow(10, v);
#endif
					else
					{
						if (p < 4)
							br.params[absValue][p] = pow(10, v);
						else
							br.params[absValue][p] = v;
					}
				}
			}
#endif
			//	Find the absolute distance from the selected 'result' LL in order
			orderedVec distanceFromOptimalLL;
			for (auto i = orderedResults.begin(); i != orderedResults.end(); i++)
				distanceFromOptimalLL.add(abs(i->first - br.LLresult));

			//	and then find the median = MAD
			br.LL_dev = distanceFromOptimalLL.median();

			//	Now go through the parameters creating ordered lists of the absolute distance from the selected param
			//	and also the positive and negative distances so that we can do single sided deviation measures
			typedef vector<orderedVec> diffList;
			vector<diffList> diffs(4, diffList(Nparams));

			size_t n = 0;
			auto j = orderedResults.begin();
			for (; n < orderedResults.size(); j++, n++)
			{
				for (size_t p = 0; p < Nparams; p++)
				{
					VEC_DATA_TYPE v = (*j->second)[p];
					if (v > br.params[logValue][p])
					{
						diffs[minLog][p].add(v - br.params[logValue][p]);
						diffs[logValue][p].add(v - br.params[logValue][p]);
					}
					else
					{
						diffs[maxLog][p].add(br.params[logValue][p] - v);
						diffs[logValue][p].add(br.params[logValue][p] - v);
					}
					if (p == 1)
#ifdef ABS_H_PARAM
						diffs[absValue][p].add(abs(br.params[absValue][p] - v));
#else
						diffs[absValue][p].add(abs(br.params[absValue][p] - pow(10, v)));
#endif
					else
					{
						if (p < 4)
							diffs[absValue][p].add(abs(br.params[absValue][p] - pow(10, v)));
						else
							diffs[absValue][p].add(abs(br.params[absValue][p] - v));
					}
				}
			}
			//	And then find the medians
			for (size_t p = 0; p < Nparams; p++)
			{
				for (size_t t = 0; t < diffs.size(); t++)
					br.param_dev[t][p] = diffs[t][p].median();
			}
			progMessage(m);
		}
	}
	/*	******************************************************************************************
	//	And then find the standard deviation
	for (modelType m :allModels())
	{
	bestResult & br = bestResults[m];
	VEC_DATA_TYPE LL_dev = 0;
	size_t size = br.params.size();
	vector<dataVec> param_dev(2, dataVec(size));
	vector<dataVec> p_N(2, dataVec(size));
	int N = 0;

	for (mcmcRunId i = 1; i <= fullResultSSChain[m].size(); i++)
	{
	for (size_t j = fullResultSSChain[m][i].size() - 1; j > fullResultSSChain[m][i].size() - 100 ; j--)
	{
	N++;
	VEC_DATA_TYPE diff = fullResultSSChain[m][i][j] - br.minLL;
	LL_dev += (diff * diff);
	for (size_t k = 0; k < size; k++)
	{
	VEC_DATA_TYPE diff = fullResultChain[m][i][j][k] - br.params[k];
	if (diff >= 0)
	{
	param_dev[0][k] += (diff * diff);
	p_N[0][k]++;
	}
	else
	{
	param_dev[1][k] += (diff * diff);
	p_N[1][k]++;
	}
	}
	}
	}
	br.minLL_dev = sqrt(LL_dev / N);
	br.param_dev[0] = sqrt(param_dev[0] / p_N[0]);
	br.param_dev[1] = sqrt(param_dev[1] / p_N[1]);
	}*/
#endif  // ifdef FIXED_RESULTS

	return true;
}

#ifdef SELECT_BY_LL
//	Find which model performed best based on the Log Liklihood
modelType LiBiNormCore::getBestModel()
{
	modelType bestModel = noModelSpecified;
	double bestLL = MAX_DOUBLE;
	for (modelType m : allModels())
	{
		if (bestResults[m].LLresult < bestLL)
		{
			bestModel = m;
			bestLL = bestResults[m].LLresult;
		}
	}
	return bestModel;
}
#endif

//	The top level summary of the results, showing best LL and associated paremeters for each model
void LiBiNormCore::printResults()
{
	string filename(outputFileroot.replaceSuffix("_results.txt"));
	TsvFile mcmcResult;
	//	Now output a table with the end points of each of the chains.
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	First headers up to and including the maximum model that is run.   Always leave space
	//	for the intermediate models so the layout of the results is consistent
	mcmcResult.printStart("");
	for (modelType m : allModels())
	{
		string count;
		if (bestResults[m].nelderMeadIterations)
			count = _s("NM iterations:", bestResults[m].nelderMeadIterations);
		mcmcResult.printMiddle(m, count);
		mcmcResult.printRepeat(headers[m].size());
	}
	mcmcResult.printEnd();
	mcmcResult.printStart("Name");

	for (modelType m : allModels())
		mcmcResult.printMiddle(headers[m], "LL", "");
	mcmcResult.printEnd();

	mcmcResult.printStart("Model");

	for (modelType m : allModels())
	{
		mcmcResult.printRepeat(headers[m].size() + 1, conv(m).substr(6).c_str());
		mcmcResult.printMiddle("");
	}
	mcmcResult.printEnd();

	//	A row for the optimal parameters that were found for each model
	for (size_t i = 0; i < 2; i++)
	{
		mcmcResult.printStart((i == 0) ? "Log Opt" : "Abs Opt");
		for (modelType m : allModels())
		{
			if (bestResults[m].params[i].size())
				mcmcResult.printMiddle(bestResults[m].params[i], bestResults[m].LLresult, "");
			else
				mcmcResult.printRepeat(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}

	//	And the initial Values
	mcmcResult.printStart("Initial");
	for (modelType m : allModels())
	{

		if (nelderMead && (initialValues[m].size()))
			mcmcResult.printMiddle(initialValues[m], "");
		else
			mcmcResult.printRepeat(headers[m].size() + 2);
	}

	mcmcResult.printEnd();

	//	A row for the deviations for each model
	for (size_t i = 0; i < 4; i++)
	{
		switch (i)
		{
		case 0:	mcmcResult.printStart("Spread Log"); break;
		case 1:	mcmcResult.printStart("Spread Abs"); break;
		case 2:	mcmcResult.printStart("pos Spread Log"); break;
		case 3:	mcmcResult.printStart("neg Spread Log"); break;
		}
		for (modelType m : allModels())
		{
			if (bestResults[m].param_dev[i].size())
				mcmcResult.printMiddle(bestResults[m].param_dev[i], bestResults[m].LL_dev, "");
			else
				mcmcResult.printRepeat(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}

	//	Print out the enpoints of all of the mcmcruns
	for (mcmcRunId i = 1; i <= Nruns; i++)
	{
		mcmcResult.printStart(_s("Chain end ", i));
		for (modelType m : allModels())
		{
			//	Was there a jth run of this model?  If so then print the end points.  *...rbegin() gets the last entry in the list.
			if (fullResultSSChain[m].size() && (i <= fullResultSSChain[m].rbegin()->first))
				mcmcResult.printMiddle(*fullResultChain[m][i].rbegin(), *fullResultSSChain[m][i].rbegin(), "");
			else
				mcmcResult.printRepeat(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	mcmcResult.close();
}



// #define BEST_RESULT_LOCATION  Use this to print out locations where best loglilihood is obtained in mcmc run
#ifdef BEST_RESULT_LOCATION
#define _BRL(A,B) A,B,
#else
#define _BRL(A,B)
#endif

void LiBiNormCore::printBias()
{
	TsvFile mcmcResult;
	string filename(outputFileroot.replaceSuffix("_norm.txt"));
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	//	Create a vector containing the list of frequencies that we are going to calculate the
	//	bias figures for
	dataVec lengths;
	lengths.push_back(DEFAULT_NORMALISATION_GENE_LENGTH);
	for (size_t i = 100; i < (size_t)min(400, MAX_GENE_LENGTH_FOR_NORM_PLOT); i += 20)
		lengths.push_back(i);
	for (size_t i = 400; i <= MAX_GENE_LENGTH_FOR_NORM_PLOT; i += 100)
		lengths.push_back(i);

	map<size_t, dataVec> biases;

	dataVec plottedLengths(lengths.subset(1));

	for (modelType m : allModels())
	{
		if (bestResults[m])
			getBias(m, bestResults[m].params[logValue], lengths, biases[m]);
	}

	for (modelType m : allModels())
	{
		mcmcResult.print(m, "Log likelihood", _BRL("mcmc run", "mcmc iteration")headers[m]);
		mcmcResult.print("Parameters", bestResults[m].LLresult, _BRL(bestResults[m].run, bestResults[m].pos) bestResults[m].params[logValue]);
		if (biases[m].size())
		{
			mcmcResult.print("Length", plottedLengths);
			mcmcResult.print("Bias", biases[m].subset(1));
		}
		else
		{
			mcmcResult.print("Length", plottedLengths);
			mcmcResult.print("Bias",dataVec(plottedLengths.size(),0));
		}
		mcmcResult.print();
	}
	mcmcResult.close();
}

#define BIN_COUNT 100

void LiBiNormCore::printDistribution(GeneCountData & geneCounts)
{
	TsvFile distResult;
	string filename(outputFileroot.replaceSuffix("_distribution.txt"));
	if (!distResult.open(filename))
		exitFail("Unable to open output File ", filename);

	vector<int> lengths{ 200,500,1000,2000,3000,4000,6000,8000,10000,15000,20000 };

	vector<dataVec> counts(lengths.size(), dataVec(BIN_COUNT));
	geneCounts.getDistribution(lengths, BIN_COUNT, counts);

	for (size_t i = 0; i < lengths.size(); i++)
	{
		counts[i].smooth(2);
		counts[i].normalise();
	}

	vector<int> entries = { 1,2,3,5,7 };

	distResult.print("ModelCount",Nmodels);
	distResult.printStart("Position");
	for (int i : entries)
	{	
		for (int j = 0;j < counts[i].size();j++)
			distResult.printMiddle((double)j/counts[i].size());
	}

	distResult.print();
	distResult.printStart("Reads");
	distResult.printEnd(counts[1], counts[2], counts[3], counts[5], counts[7]);
	distResult.printStart("Length");
	for (int i : entries)
		distResult.printRepeat(counts[i].size(), lengths[i]);
	distResult.printEnd();


	if (Nmodels)
	{
		distResult.printStart(theModel);

		for (int i : entries)
		{
			dataVec dist = getDistribution(theModel, lengths[i], BIN_COUNT, bestResults[theModel].params[logValue]);
			distResult.printMiddle(dist);
		}
		distResult.printEnd();

		for (modelType modl : allModels())
		{
			if (bestResults[modl].params[logValue].size())
			{
				distResult.printStart(_s(modl, " LL=", $("%7.0f", -bestResults[modl].LLresult)));
				for (int i : entries)
				{
					dataVec dist = getDistribution(modl, lengths[i], BIN_COUNT, bestResults[modl].params[logValue]);
					distResult.printMiddle(dist);
				}
			}
			else
			{
				distResult.printStart(modl);
				distResult.printRepeat(BIN_COUNT * entries.size());
			}
			distResult.printEnd();
		}
		distResult.print();
	}

	distResult.print("Reads");
	for (size_t i = 0; i < lengths.size(); i++)
			distResult.print(lengths[i], counts[i]);
	distResult.print();

	for (modelType modl : allModels())
	{
		if (bestResults[modl].params[logValue].size())
		{
			distResult.print(modl);

			for (auto l : lengths)
			{
				dataVec dist = getDistribution(modl, l, BIN_COUNT, bestResults[modl].params[logValue]);
				distResult.print(l, dist);
			}
			distResult.print();
		}
	}
}


void LiBiNormCore::printAllMcmcRunData()
{
	TsvFile mcmcResult;
	for (modelType modl : allModels())
	{
		string filename = outputFileroot.replaceSuffix("_", conv(modl, true), ".txt");
		if (!mcmcResult.open(filename))
			exitFail("Unable to open output file ", filename);

		//			mcmcResult.printMiddle(headers[modl], "chain", "");

		//	This ensures that at least one header is output, which ensures that there is something in the file
		//	even if no data were produced for this model
		for (size_t i = 0; i < fullResultChain[modl].size(); i++)
			mcmcResult.printMiddle(headers[modl], "LL", "");
		mcmcResult.printEnd();

		//	For each of the mcmc runs print out the results.  Each run is a column
		for (size_t i = 0; i < fullResultChain[modl][1].size(); i++)
		{
			for (mcmcRunId j = 1; j <= fullResultChain[modl].rbegin()->first; j++)
			{
				// fullResultChain[modl][j][i] is a vector of N values which are printed out as N tab separated values
				//	using the TsvFile support for printing vectors.
				mcmcResult.printMiddle(fullResultChain[modl][j][i], fullResultSSChain[modl][j][i], "");
			}
			mcmcResult.printEnd();
		}
		mcmcResult.close();
	}
}

void LiBiNormCore::printConsolidatedMcmcRunData()
{
	TsvFile mcmcResult;
	string filename = outputFileroot.replaceSuffix("_model_cons.txt");
	if (!mcmcResult.open(filename))
		exitFail("Unable to open output File ", filename);

	map<size_t, multimap <double, dataVec *>::iterator > iterators;
	mcmcResult.printStart("");
	for (modelType m : allModels())
	{
		mcmcResult.printMiddle(headers[m], "chain", "");
		iterators[m] = allOrderedResults[m].begin();
	}
	mcmcResult.printEnd("");
	bool found = true;
	for (size_t i = 0; (i < 1000000) && found; i++)
	{
		found = false;
		mcmcResult.printStart(i);
		for (modelType m : allModels())
		{
			if (iterators[m] != allOrderedResults[m].end())
			{
				found = true;
				mcmcResult.printMiddle(*(iterators[m]->second), iterators[m]->first, "");
				iterators[m]++;
			}
			else
				mcmcResult.printRepeat(headers[m].size() + 2);
		}
		mcmcResult.printEnd();
	}
	mcmcResult.close();
}
