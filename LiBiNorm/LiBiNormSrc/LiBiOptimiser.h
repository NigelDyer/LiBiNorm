// ***************************************************************************
// LiBiOptimiser.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Extends the core NelderMead optimiser for paremeter estimation 
// ***************************************************************************

#ifndef LIBIOPTIMISER_H
#define LIBIOPTIMISER_H

#include "nelderMeadOptimiser.h"
#include "ModelData.h"
#include "mcmc.h"

class LiBiOptimiser : public nelderMeadOptimiser
{
public:
	LiBiOptimiser(const mcmcGeneData & geneData, modelType currentModel,int & iterations) : 
		nelderMeadOptimiser(iterations),
		geneData(geneData), currentModel(currentModel) {};
	virtual ~LiBiOptimiser(){};
public:


	dataVec getParams(modelType m,optionsType & options,dataVec & initialValues);

	//	Generates the error/liklyhood/sum of squares for the current model + paraneter values
	virtual ErrorPair ErrorFunc();
	//	Can be used to print results to screen/file during the process.  Currently ignored
	virtual void SaveResults(bool toFile) {};

	int iterations;

	//	
	optionsType * opts;
	const mcmcGeneData & geneData;
	const modelType currentModel;

private:
	optiDataType allOptiData;
	paramDescriptionSet params;

};

#endif

