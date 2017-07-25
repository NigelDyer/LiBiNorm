// ***************************************************************************
// rand.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some packaged C++11 code for generating random numbers
// ***************************************************************************

#include "rand.h"
#include <random>


//  Returns a random number, uniformly distributed between 0 and a.
double rand(double a)
{
	// Static members generate distribution in range 0 to 1.   This is scaled by input parameter for each call
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<double> dis(0, 1);

	return dis(gen) * a;
}

//	Sets up a vector of values with a normal distribution
dataVec randn(size_t x)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static	std::normal_distribution<double> distribution;

	dataVec retVal(x);
	for (auto & i : retVal)
		i = distribution(gen);
	return retVal;
}
