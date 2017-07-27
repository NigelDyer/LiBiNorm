// ***************************************************************************
// rand.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some packaged C++11 code for generating random numbers
// ***************************************************************************

#ifndef RAND_H
#define RAND_H

#include "dataVec.h"

//  Returns a random number, uniformly distributed between 0 and a.
double rand(double a);

//	Sets up a vector of values with a normal distribution
dataVec randn(size_t x);

#endif