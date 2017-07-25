// ***************************************************************************
// ModelData.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Code for each of the five models
// ***************************************************************************

#ifndef MODELDATA_H
#define MODELDATA_H

#include <float.h>
#include <map>
#include <algorithm>
#include <vector>
#include "mcmc.h"

/*
Holds the information relating to the 6 different models 
*/

//	The enumerated types used globally to identify the models
enum modelType
{
	ModelA = 0,
	ModelB = 1,
	ModelC = 2,
	ModelD = 3,
	ModelE = 4,
	ModelBD = 5,
	noModel = 6,
	findBestModel = 7
};

//	A collection of all the models, allowing code to iterate through them
const std::vector<modelType> & allModels();

void setSSfun(optionsType & options, modelType m);

//	Routines for converting between model identifier and equivalent strings
std::string conv(const modelType m,bool removeGaps = false);
modelType modelFromString(const std::string & desc);
//	Allows the model to be output to a stream such as std::out as appropriate text
bool printVal(outputDataFile * f, modelType m);
//	Allows genomicPositions to be placed into StringEx's
namespace std
{
	std::string to_string(const modelType & m);
}
std::ostream& operator<< (std::ostream &out, const modelType & m);

inline void parseval(const char * start, modelType & value, size_t & len)
{
	value = modelFromString(std::string(start, len));
};



//  Returns mcmc paremetyers associated with a model
paramSet GetModelParams(modelType model,dataVec * defaults = 0,VEC_DATA_TYPE offset = 0);

//	The set of headers associated with the model parameters, extracted from the data provided 
//	by GetModelParams
typedef std::map<modelType, std::vector<std::string> > headerType;
headerType getHeaders();


//	Holds the information relating to the best set of parameters associated with each model
struct bestResult
{
	bestResult() :LLresult(DBL_MAX), LL_dev(0), run(0), pos(0), nelderMeadIterations(0) {};
	//	For checking if paraeters have been estimated
	operator bool() const { return params[0].size(); };
	// The miniumum log liklyhood, and the associated deviation
	VEC_DATA_TYPE LLresult, LL_dev;
	//	The mcmc run, and the position in the run where the best result was found
	size_t run, pos;

	int nelderMeadIterations;
	//	The best model parameters found
	dataVec params[2];
	//	A measure of the positive and negative spread of the parameters and the absolute diff
	dataVec param_dev[4];
};

//	Gets the bias for a selection of lengths
void getBias(modelType m, dataVec & params, const dataVec & lengths,dataVec & bias);

//	Calulates the log liklyhoods for each model
double FLL_ModelA(const dataVec & param, const mcmcGeneData & data);
double FLL_ModelB(const dataVec & param, const mcmcGeneData & data);
double FLL_ModelC(const dataVec & param, const mcmcGeneData & data);
double FLL_ModelD(const dataVec & param, const mcmcGeneData & data);
double FLL_ModelE(const dataVec & param, const mcmcGeneData & data);
double FLL_ModelBD(const dataVec & param, const mcmcGeneData & data);


// Common prior function
double priorFunc(const dataVec & data, const paramSet & params);
double priorFuncE(const dataVec & data, const paramSet & params);
#endif

