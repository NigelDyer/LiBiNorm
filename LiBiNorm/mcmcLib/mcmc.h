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
#include "nelderMeadOptimiser.h"

/*
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

	double(*ssfun)(const dataVec & param, const mcmcGeneData & data);
	double(*priorfun)(const dataVec &, const paramSet &);
	double sigma2;
	size_t nsimu,Nruns;
	double jumpSize;
	dataVec qcov;
};


//	The class that manages the Monte Carlo Markov Chain determintaion of parameters
class mcmc
{
public:
	void mcmcrun(const mcmcGeneData & data,const paramSet & params,const optionsType & options);
	
	double priorfun(const dataVec & th, const dataVec & mu, const dataVec & sig);

	std::vector<dataVec> _chain;
	dataVec _sschain;
};


#endif // !MCMC_H

