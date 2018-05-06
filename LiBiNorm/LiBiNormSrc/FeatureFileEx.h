// ***************************************************************************
// FeatureFileEx.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 21 April 2018
// ---------------------------------------------------------------------------
// Extends featureFile class in bioinformaticsLib to support processing of 
// gff and gtf files
// ***************************************************************************

#ifndef FEATURE_FILE_EX_HEADER
#define FEATURE_FILE_EX_HEADER

#include "containerEx.h"
#include "libCommon.h"
#include "Regions.h"
#include "featureFile.h"
#include "dataVec.h"
#include "GeneCountData.h"

struct featureOverlap
{
	featureOverlap(size_t start, size_t finish, rna_pos_type RNAstartPos, rna_pos_type RNAendPos, bool strict, const std::string & geneName, const std::string & featType) :
		start(start), finish(finish), RNAstartPos(RNAstartPos), RNAendPos(RNAendPos), strict(strict), geneName(geneName), featType(featType)
	{};

	size_t start,finish;
	rna_pos_type RNAstartPos,RNAendPos;
	bool strict;
	//	Store references here as these are created and deleted lots and these removes the need to allocate and
	//	deallocate on the heap.
	const std::string & geneName;
	const std::string & featType;
};

class featureRegion
{
public:
	class  chromosomeFeatureData : public std::multimap<size_t, featureRegion>
	{
	public:
		ADD_ITER(position, feature);
		ADD_CONSTMAPPAIR_VAR(position, feature);
	};

	featureRegion(featureRegion && gtf);
	featureRegion(size_t start, size_t finish, const std::string & name, char strand, const std::string & type, const std::string & bioType);
	
	void checkOverlap(const region & segment, std::vector<featureOverlap> & overlapList) const;

	size_t start,finish;
	rna_pos_type RNAstart;
	//	Store actual values here as these are only created once and then referenced lots, so this is more efficient
	const std::string name;
	const std::string type;
	const std::string bioType;
	char strand;
	chromosomeFeatureData::iterator overlaps;

#ifdef HISAT2
	std::string sequence;
#endif
};


//	A map of the features associated with each gene
class genomeFeatureRegions : public std::map<std::string, featureRegion::chromosomeFeatureData>
{

};

//	For each region in the genome where a read starts we store a reference to an iterator that 
//	points to the first region that overlaps the region 
class chromosomePositionIndexMap : public std::multimap<size_t, featureRegion::chromosomeFeatureData::Iterator>
{
public:
	ADD_ITER(position,featureDataIterator)
};

typedef chromosomePositionIndexMap chromosomeEndIndexMap;
typedef std::map<std::string,chromosomeEndIndexMap > genomeEndIndexMap;

typedef std::vector<featureRegion *> featureRegionList;

class geneData
{
public:
	geneData() : overlapsAnotherGene(false), strand(' ') {};

	void addRegion(featureRegion * newRegion,bool ol, VEC_DATA_TYPE & length,const std::string & chr)
	{
		if (ol)
			overlapsAnotherGene = true;

		newRegion->RNAstart = length + 1;
		if (regions.size() == 0)
		{
			strand = newRegion->strand;
			chromosome = chr;
		}

		regions.push_back(newRegion);
		length += (newRegion->finish - newRegion->start + 1);

	};

	featureRegionList regions;
	std::string chromosome;
	bool overlapsAnotherGene;
	char strand;

#ifdef HISAT2
	size_t length;
	std::string priorSeq, postSeq;
	std::map<rna_pos_type, size_t> mRNAtoSeq;
#endif

};

//
//	The featureFile class does the basic parsing of a gtf or gff class.  featureFileEx extends this to 
//	provide the specific information required for identifying reads with genes in LiBiNorm count
class featureFileEx : public featureFile
{
public: 
	//  Create an index of the features once they have been read from a file
	void index(GeneCountData & geneCounts,bool useStrand);
	//	Print feature data and counts to a file
	bool outputChromData(const std::string & filename, const GeneCountData & geneCounts);
	//	Print feature data as bed file
	bool outputBedData(const std::string & filename);

	//	A container of all the consolidated feature regions
	genomeFeatureRegions genomeGtfData; 
	//	A map of the ends of the feature regions.   Used for finding overlaps
	genomeEndIndexMap genomeEndIndex;
	//	A map of the regions associated with a gene, indexed by the unique identifier of the gene (which could be a transcript)
	std::map<std::string,geneData> genes;
};

#endif
