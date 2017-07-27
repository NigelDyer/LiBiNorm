// ***************************************************************************
// mcmc.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// An implementation of the mcmc algorithm
// ***************************************************************************

#ifndef MCMC_H
#define MCMC_H

#include "params.h"

/*
	Set of options for parameter estimation.   The ssfun and priorFun methods are used both 
	for the initial Nelder Mead stage of the optimisation and also mcmc stage

% options structure
%    options.ssfun    -2*log(likelihood) function
%    options.sigma2   initial error variance
%    options.Nsimu    number of mcmc cycles in a run
%    options.Nruns    total number of runs/observations
*/
class optionsType
{
public:	
	optionsType():nsimu(3){};

	//	Sum of squares error function given a set of parameters
	double(*ssfun)(const dataVec & param, const mcmcGeneData & data);
	//	Prior given a set of parameters
	double(*priorfun)(const dataVec &, const paramDescriptionSet &);

	//	mcmc related values
	double sigma2;
	size_t nsimu,Nruns;
	double jumpSize;
	dataVec qcov;
};


//	The class that manages the Monte Carlo Markov Chain determintaion of parameters
class mcmc
{
public:
	void mcmcrun(const mcmcGeneData & data,const paramDescriptionSet & params,const optionsType & options);

	//	The chain of parameter values generated during the mcmc run
	std::vector<dataVec> _chain;
	//	And the associated sum of squares errors
	dataVec _sschain;
};


#endif // !MCMC_H

