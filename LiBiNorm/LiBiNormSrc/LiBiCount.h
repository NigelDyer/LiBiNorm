// ***************************************************************************
// LiBiCount.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// The top level code associated with "LiBiNorm count" modes
// ***************************************************************************

#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "FeatureFileEx.h"
#include "api/BamWriter.h"
#include "LiBiNorm.h"

using namespace BamTools;


//	The htseq-count modes
enum mode {
	mode_none,
	intersect_union,
	intersect_strict,
	intersect_nonempty,
	intersect_all
};



class LiBiCount : private LiBiNormCore
{
	//	Used for reading back cached read information from cache files
	class cacheEntry
	{
	public:
		cacheEntry() : file(0) {};
		~cacheEntry();

		bool open(const std::string filename);
		bool readNext();
		void close();

		//	Holds the name of the read, which is not in the readData class
		readData currentRead;
		std::string name;
	private:
		std::string fname;
		std::ifstream * file;
	};



public:
	LiBiCount() : bamOutMode(outputNone), landscapeFile(false), countMode(mode_none) {};

	int main(int argc, char **argv);

private:
	//	For reading and processing bam data
	bool processNameOrderedBamData();
	bool processPositionOrderedBamData();
	void processCachedReads(size_t cacheFileCount);

	//Support functions for reading bam data
	bool AReadIsMapped(const BamAlignment & ba);
	std::string addRead(const regionLists & segments, const featureFileEx & gtfData);
	void incBamCounter(const BamAlignment * ba = 0, int size = -1);

	//	For comparing two files.  Not currently used
	void fileCompare(int argc, char **argv);
	
	//  For reading and persisting header data from a bam file
	BamReader reader;
	BamWriter writer;
	enum BamOutMode
	{
		outputNone,
		outputAll,
		outputMatched,
		outputUnmatched
	} bamOutMode;
	RefVector references;

	//	Config data for reading the bam file
	bool useStrand, reverseStrand, nameOrder,landscapeFile;
	mode countMode;
	int minqual;
	size_t bamCounter, maxCacheSize;

	//	Directory data for outputting results
	stringEx tempDirectory;
	TsvFile genomeDataFile;

	//	For obtaining information from the gtf/gff file
	featureFileEx genomeDef;
};


#endif
