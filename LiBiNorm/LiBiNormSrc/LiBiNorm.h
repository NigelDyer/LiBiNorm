// ***************************************************************************
// LiBiNorm.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// The top level code associated with "LiBiNorm model" mode
// ***************************************************************************

#ifndef LIBINORM_H
#define LIBINORM_H

#include "LiBiNormCore.h"

class LiBiNorm : protected LiBiNormCore
{
public:
	int main(int argc, char **argv);

private:
	stringEx landscapeFilename;
};

#endif


