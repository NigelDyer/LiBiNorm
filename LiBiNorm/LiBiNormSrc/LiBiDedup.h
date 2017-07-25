// ***************************************************************************
// LiBeDedup.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For removing duplicate entries in BamFiles
// ***************************************************************************

#pragma once
class LiBiDedup
{
	int minCount,windowSize;

public:
	int main(int argc, char **argv);

};

