// ***************************************************************************
// Regions.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For information relating to regions
// ***************************************************************************

#include "Regions.h"

using namespace std;
using namespace BamTools;


//	This constructor creates the readData from the bam file entry.  This means that methods expecting 
//	readData can be passed a bamAlignment.  Use an rValue constructor so that we can 'swallow up' the cigar data 
//	rather than making a copy of it as once the readData has been created we will have no further use
//	for the cigar data
readData::readData(BamTools::BamAlignment && ba) :
	refId(ba.RefID),
	position(ba.Position + 1),
	qual(ba.MapQuality),
	cigar(move(ba.CigarData))
{
	//	This simulates line 155 in count py.  If there is no optional NH field in the first read
	//	then the python code throws an error, which we simulate by setting NH to zero
	if (!ba.GetTag("NH", (unsigned short &)NH))
		NH = ba.IsFirstMate() ? 0 : -1;
	if (ba.IsPaired())
	{
		if (ba.IsReverseStrand() == ba.IsFirstMate())
			strand = '-';
		else
			strand = '+';
	}
	else
		strand = ba.IsReverseStrand() ? '-' : '+';
};

readData::readData(const BamTools::BamAlignment & ba) :
	refId(ba.RefID),
	position(ba.Position + 1),
	qual(ba.MapQuality),
	cigar(ba.CigarData)
{
	//	This simulates line 155 in count py.  If there is no optional NH field in the first read
	//	then the python code throws an error, which we simulate by setting NH to zero
	if (!ba.GetTag("NH", (unsigned short&)NH))
		NH = ba.IsFirstMate() ? 0 : -1;
	if (ba.IsPaired())
		strand = (ba.IsReverseStrand() == ba.IsFirstMate()) ? '-' : '+';
	else
		strand = ba.IsReverseStrand() ? '-' : '+';
};


//	Combines a region with an existing set of regions
void regionList::combineRegion(const region & r1)
{
	bool combined = false;

	for (auto i = data.begin();i != data.end();i++)
	{
		//	Does the new region start before the end of the existinng region, and is it on the same strand?
		if ((r1.start <= i->second.end) && (r1.strand == i->second.strand))
		{
			//	If it starts before the existing start and ends after it then there is an overlap 
			//	and we are going to have to replace the existing region as regions are indexed by
			//  the start
			if ((r1.start <= i->second.start) && (r1.end >= i->second.start)) 
			{
				// The new region, based on the existing one
				region r = i->second;
				// but with an earlier start
				r.start = r1.start;
				//	and possibly an earlier end.
				if (r1.end >= i->second.end)
					r.end = r1.end;

				//	replace the existing entry with the new one, taking care not to saw off the
				//	branch you are sitting on and 
				data.erase(i);
				i = data.emplace(r1.start,r);
				combined = true;
				break;	
			}
			//	The new region just extends the end of the existing region
			else if (r1.end >= i->second.end)
			{
				i->second.end = r1.end;
				combined = true;
				break;	
			}
		}
	}
	//	If this is a new region then add it to the map
	if (!combined)
		data.emplace(r1.start,r1);
}

//	Combine this region list with another (they are on the same chromosome
void regionList::combine(const regionList & rl)
{
	for (auto r : rl.data)
		combineRegion(r.second);
}

//	Constructs a set of regions from the information in a read, ie bam object
regionList::regionList(const readData & read)
{
	// initialize alignment end to starting position
	size_t start = read.position;
	size_t end = start;

	// iterate over cigar operations
	for (const CigarOp & co : read.cigar) {

		switch ( co.Type ) {

			// increase end position on CIGAR chars [DMXN=]
			case Constants::BAM_CIGAR_DEL_CHAR      :
			case Constants::BAM_CIGAR_MATCH_CHAR    :
			case Constants::BAM_CIGAR_MISMATCH_CHAR :
			case Constants::BAM_CIGAR_SEQMATCH_CHAR :
				end += co.Length;
				break;

			case Constants::BAM_CIGAR_INS_CHAR :
				break;

			//	If there is a skip then the region after the skip is a new region
			case Constants::BAM_CIGAR_REFSKIP_CHAR  :
				{
					combineRegion(region(start,end-1,read.strand));
					start = (end + co.Length);
					end = start;
					break;
				}
			default :
				break;
		}
	}
	combineRegion(region(start,end-1,read.strand));
}

//	Creates a regionList from one of the reads, either from a bam entry or from cachedData.  Use emplace so that the
//	regionList can be efficiently placed straight into the map.
regionLists::regionLists(const readData & read, const std::string & name) :
	name(name), NH(read.NH),
	qual(read.qual),
	strands(1, read.strand)
{
	data.emplace(read.refId, regionList(read));
};

//	Adds the information associated with the second read, which will be placed in the existing chromosome
//  or added to a new.
void regionLists::combine(const readData & read) 
{
	data[read.refId].combine(regionList(read));
	if ((read.NH == 0) || (NH == 0))
		NH = 0;
	else if (read.NH > NH)
		NH = read.NH;
	if (read.qual < qual)
		qual = read.qual;
	strands.emplace_back(read.strand);
};




//	Parse a text string and convert it into a cigar value.  Used when retrieving entries from cache files
void parserInternal::parseval(const char *& start,Cigar & co,size_t & len)
{
	size_t i = 0;
	while (i < len)
	{
		char c = start[i++];
		int val = 0;
		while ((start[i] >= '0') && (start[i] <= '9') && (i < len))
			val = (val *10)+ (start[i++]-'0');
		co.emplace_back(CigarOp(c,val));
	}
}

//	Prints out cigar information, used by the printTsv class
bool printVal(outputDataFile * f,const Cigar & cigar)
{
	for (auto & i: cigar)
		fprintf(f->fout,"%c%i",i.Type,i.Length);
	return true;
};

//	Prints the value of a readData instance.  Note that printing read.cigar will have the effect of calling the 
//	printVal associated with the cigar above
bool printVal(outputDataFile * f,const readData & read)
{
	f->printStart(read.refId,read.position,read.strand,read.cigar,read.NH,read.qual);
	return true;
};


