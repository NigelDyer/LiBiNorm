// ***************************************************************************
// LiBiOptimiser.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Extends the core NelderMead optimiser for paremeter estimation 
// ***************************************************************************
#include "rand.h"
#include "Options.h"
#include "LiBiOptimiser.h"

using namespace std;

ErrorPair LiBiOptimiser::ErrorFunc()
{
	ErrorPair ep;
	dataVec data;
	for (size_t i = 0; i < (*allOptiData[0]).size(); i++)
		data.push_back((*allOptiData[0])[i].value);

	ep = opts->ssfun(data, geneData);
	ep.WithPrior += opts->priorfun(data, params);

/*	cout << ep.WithPrior << endl;
	string s;
	cin >> s;*/
	return ep;
};


dataVec LiBiOptimiser::getParams(modelType m, optionsType & options, dataVec & initialValues)
{
	opts = & options;
	params  = GetModelParams(m);

	optiVector optiData;
	for (size_t i = 0;i < params.size();i++)
	{
		if (initialValues.size() > i)
			optiData.push_back(optiItem(initialValues[i], true));
		else
			optiData.push_back(optiItem(params[i].initial, true));
	}

	allOptiData.push_back(&optiData);
	if (params.size() > 2)
	{
		for (size_t i = 0; i < 2; i++)
			allOptiData[0][0][i].optimise = false;

		optimise(allOptiData, 100, 20, options.jumpSize, conv(m));
		optimise(allOptiData, 100, 20, options.jumpSize, conv(m));

		for (size_t i = 0; i < 2; i++)
			allOptiData[0][0][i].optimise = true;

		optimise(allOptiData, 100, 20, options.jumpSize, conv(m));
	}

	VEC_DATA_TYPE LL = optimise(allOptiData, NELDER_MEAD_ITERATIONS,20, options.jumpSize,conv(m));

	dataVec results;
	for (size_t i = 0; i < params.size(); i++)
		results.push_back(allOptiData[0]->at(i).value);
	results.push_back(LL);
	return results;
}
