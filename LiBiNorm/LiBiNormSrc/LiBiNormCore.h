// ***************************************************************************
// LiBiNormCore.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// Common code associated with all of the LiBiNorm modes
// ***************************************************************************

#ifndef LIBINORMCORE_H
#define LIBINORMCORE_H

#include <map>
#include <mutex>
#include "stringEx.h"
#include "containerEx.h"
#include "GeneCountData.h"
#include "Options.h"
#include "ModelData.h"

//	For identifying the mcmc run for a specific model.  Numbered from 1.
typedef size_t mcmcRunId;


//	LiBiNormCore contains common functionality for LiBiNorm count. model and variation
class LiBiNormCore
{
public:
	//	Used to identify param and param_dev values within bestResults
	enum
	{
		logValue = 0,
		absValue = 1,
		minLog = 2,
		maxLog = 3
	};


protected:
	LiBiNormCore() :normalise(true), pauseAtEnd(false), 
#ifdef USE_NELDER_MEAD_FOR_INITIAL_VALUES
		nelderMead(true),
		Nsimu(NELDER_MCMC_ITERATIONS),
#else
		nelderMead(false),
		Nsimu(MCMC_ITERATIONS),
#endif
#ifdef SELECT_BY_LL
		theModel(noModelSpecified), bestModel(noModelSpecified),
#else
		theModel(ModelBD), calcAllModels(false),
#endif
		maxReads(DEF_MAX_READS_FOR_PARAM_ESTIMATION),
		maxGeneLength(DEF_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION),
		Nthreads(DEF_THREADS),
		Nruns(NUMBER_OF_MCMC_RUNS),
		NrunsOtherModels(0),
		outputFPKM(false)
	{
		headers = getHeaders();
	};

	void helpCommon();
	void featureAndIdAttributeHelp();

	bool commandParseCommon(int & ni, char **argv);
	bool commandParseIdAndType(int & ni, char **argv);

	void SetInitialParamsFromFile(const std::string & filename);
	void checkFeatureAndIdAttribute();

	void mcmcThread(optionsType options);
	bool coreParameterEstimation();

	void printResults();
	void printBias();
	void printDistribution(GeneCountData & geneCounts);
	void printAllMcmcRunData();
	void printConsolidatedMcmcRunData();


protected:

	bool normalise, pauseAtEnd,nelderMead;
	size_t Nsimu,Nmodels;
	modelType theModel;

#ifdef SELECT_BY_LL
	modelType bestModel;
	modelType getBestModel();
#else
	bool calcAllModels;
#endif

	//	Information relating to the gtf/gff file being used
	//	use by LiBiNorm count and some tools
	stringEx featureFileName;
	stringEx id_attribute, feature_type;

	size_t maxReads,maxGeneLength,Nthreads;
	mcmcRunId Nruns,NrunsOtherModels;

	headerType headers;
	stringEx outputFileroot, countsFilename, parameterFilename;

	std::map<modelType, dataVec> initialValues;


	//	These vector holds the full results for each of the models, which are needed for identifying
	//	the optimal parameters and the variation that is seen.
	//	For each model the results for each mcmc run is stored as a map indexd by run number
	//	This is because the runs are done on separate threads and we want to store the results by the run
	//	identifier and not the order that they finished
	std::map<modelType, std::map<mcmcRunId, std::vector <dataVec > > >fullResultChain;
	std::map<modelType, std::map<mcmcRunId, dataVec> >fullResultSSChain;

	//	Holds the results from multiple chains combined into a single ordered list
	std::map<modelType, std::multimap <double, dataVec *> > allOrderedResults;

	std::map<modelType, bestResult> bestResults;

	//	The specific data that will be used for the mcmc parameter determination
	mcmcGeneData allGeneData;
	GeneCountData allGeneCounts;
private:
	//	Counts of the number of mcmc runs that will be done for each model
	struct loop_counts { mcmcRunId requested, counter; };
	std::map<modelType, loop_counts> threadLoopCounts;
protected:
	bool outputFPKM;
	size_t nelderMeadCounter;

};

typedef std::map<modelType, std::vector<std::string> > headerType;


#endif


