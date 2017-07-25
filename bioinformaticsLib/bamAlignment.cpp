// ***************************************************************************
// BamAlignmentEx.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Extends bamtools to provide a sort function.  BAM_LIBRARY must be defined for this to be 
//	available
// ***************************************************************************

#ifdef BAM_LIBRARY

#include "api/BamReader.h"
#include "toolkit/bamtools_sort.h"

#include "libCommon.h" //	For strdup redefinition

#include "stringEx.h"
using namespace BamTools;
using namespace std;

bool sortBamFile(const string & filename, const std::string suffix)
{
#define ARGNO 6
	SortTool sorter;

	char * args[ARGNO];
	args[0] = strdup("");
	args[1] = strdup("");
	args[2] = strdup("-in");
	args[3] = strdup(filename.c_str());
	args[4] = strdup("-out");
	args[5] = strdup(stringEx(filename).replaceSuffix(suffix,".bam").c_str());

	sorter.Run(ARGNO,args);
	for (int i = 0;i < ARGNO;i++)
		free (args[i]);

	BamReader newBam;
	newBam.Open(stringEx(filename).replaceSuffix(suffix,".bam"));
	newBam.CreateIndex();
	newBam.Close();
	return true;
}

#endif