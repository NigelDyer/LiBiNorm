// ***************************************************************************
// Regions.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For information relating to regions
// ***************************************************************************

#ifndef REGIONS_H
#define REGIONS_H

#include<map>
#include<string>
#include "printEx.h"
#include "api/BamReader.h"
#include "Options.h"

//	Wrap the cigar data so that the parser and printVal methods will recognise it
class Cigar : public std::vector<BamTools::CigarOp>
{
public:
	Cigar(){};
	Cigar(const std::vector<BamTools::CigarOp> & co): std::vector<BamTools::CigarOp>(co) {};
};

//	This holds the data from a bam entry that we are actually interested in.  It is the basis
//	of 'in program' persisted data and also data that is persisted to cache files.  Does not 
//	include the read name as this is stored eleswhere (e.g as the index of a map of readData)
class readData
{
	public:
		readData(void):NH(0),qual(0),refId(0),position(0) 
		{}
		readData(BamTools::BamAlignment && ba);
		readData(const BamTools::BamAlignment & ba);

		int refId, position;
		char strand;
		short NH, qual;
		Cigar cigar;
};

//	A reagion within a chromosome
class region 
{
public:
	size_t start,end;
	char strand;
	region(size_t start, size_t end,char strand) : start(start),end(end),strand(strand) {};
};


//	A set of regions on one chromosome
class regionList 
{
public:
	regionList(){};

	//	Which is normally created from a read
	regionList(const readData & read);

	//	Map of the regions, indexed by the location on teh chromosome
	//	Separate entries for and -ve strands so needs to be a multimap to cater 
	//	for + and - entries starting at the same location (usually an artefact)
	std::multimap<size_t,region> data;


	//	For combining data from a second read
	void combine(const regionList & rl);

private:
		//	For combining each individual read within a regionList
	void combineRegion(const region & r);

};

//	The list of regions associated with a read pair
class regionLists
{
public:
	//	Contains the regions themselves in a map indexed by chromosome (indicated by refId, the chromosome identifier
	//	in the bam file.  Done this way to cater for a read pair where the reads are on difference chromosomes
	std::map<int,regionList> data;
	//	The name of the read
	const std::string & name;
	int NH;
	int qual;
	std::vector<char> strands;

	//	Creates a regionList from one of the reads, either from a bam entry or from cachedData.  
	regionLists(const readData & read, const std::string & name);
	
	//	Adds the information associated with the second read
	void combine(const readData & read);
	
};

//	Declare the availability of methods that are used by parser for parseing cigar strings
//	and also the methods used for printing a cacheentry and the cigar data
//	These are both used for the temporary cache data that is placed on disk
namespace parserInternal
{
	void parseval(const char *& start, Cigar & cigar,size_t & len);
}

bool printVal(outputDataFile * f,const Cigar & cigar);
bool printVal(outputDataFile * f,const readData & read);

#endif

