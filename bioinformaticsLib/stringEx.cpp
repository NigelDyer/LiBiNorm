// ***************************************************************************
// stringEx.cpp (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// Extension to std::string class
// ***************************************************************************
#include <algorithm>
#include "stringEx.h"

using namespace std;

stringEx stringEx::toLower() const {
	string retVal(*this);
	transform(retVal.begin(), retVal.end(), retVal.begin(), ::tolower);
	return retVal;
}
