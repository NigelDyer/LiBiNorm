// ***************************************************************************
// LiBiTools.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// Assorted tools associated with LiBiNorm
// ***************************************************************************

#ifndef LIBITOOLS_H
#define LIBITOOLS_H

#include "LiBiNormCore.h"

class LiBiTools : public LiBiNormCore
{
public:
	int landMain(int argc, char **argv);
	int landMain2(int argc, char **argv);
	int geneMain(int argc, char **argv);
	int refSeqsMain(int argc, char **argv);
	int bedMain(int argc, char **argv);
};

#endif // !LIBITOOLS_H


