// ***************************************************************************
// BamAlignmentEx.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Extends bamtools to provide a sort function.  BAM_LIBRARY must be defined for this to be 
//	available
// ***************************************************************************


#ifdef BAM_LIBRARY

#ifndef BAMALIGNMENTEX_H
#define BAMALIGNMENTEX_H

#include <string>
bool sortBamFile(const std::string & filename,const std::string suffix = ".sort");


#endif
#endif