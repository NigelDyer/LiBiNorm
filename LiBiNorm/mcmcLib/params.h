// ***************************************************************************
// params.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Parameter information for use by the mcmc and nelder mead parameter esztimation
//  algorithms
// ***************************************************************************

#ifndef PARAMS_H
#define PARAMS_H

#include <limits> 
#include "nelderMeadOptimiser.h"
#include "dataVec.h"


//	Holds the subset of the read information associated with a gene that is used by the functions for 
//	calculating log liklyhoods with the mcmc
class mcmcGeneData
{
public:
	//	A list of all the read positions
	dataVec fragData;
	//	The list of the indexes to the genes associated with the fragments
	std::vector<int> geneIndex;
	//	The data associated with each gene, namely the length and the frequency with which it occurs
	dataVec geneLengths;
	dataVec geneFrequencies;
};


class paramType
{
public:
	paramType(std::string name, double min = MIN_DOUBLE, double max = MAX_DOUBLE, double initial = 0);
	std::string name;
	double min, max, value,initial;
	bool targetflag, localflag;
};

class paramSet :public std::vector<paramType>
{
public:
	paramSet() {};
	//	This = operator allows the values to be assigned using an aggregate initiailiser for the base class
	paramSet & operator = (std::vector<paramType>(a)) { std::vector<paramType>::operator = (a); return *this; };

	bool isValid(const dataVec & data) const;
	void setValues(const dataVec & vals);

	dataVec getvalues() const;

};

#endif