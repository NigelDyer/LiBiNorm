// ***************************************************************************
// GeneCountData.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 30 March 2018
// ---------------------------------------------------------------------------
// For processing count information associated with genes
// ***************************************************************************

#ifndef GENE_COUNT_DATA_H
#define GENE_COUNT_DATA_H

#include <random>
#include "containerEx.h"

#include "libParser.h"
#include "mcmc.h"


//
//	Used to generate random integers between 0 and max.  This does not use the C++ 
//	function to generate a flat distribution because its implementation (unlike the
//	underlying random number generator) is not consistent between compilers.  This consistency
//	is needed for comparison testing
class intRandClass
{
private:
	intRandClass()
	{
#ifdef SELECT_READS_SEED
		seed = SELECT_READS_SEED;
#else
		seed = rd();
#endif
		gen.seed(seed);
	};

public:
	//	A single instance is used for generating random numbers
	static intRandClass & instance();
	//	Use simple modulus function to generate a number between 0 and max-1
	unsigned int value(unsigned int max);
	//	Set a specific seed
	void reseed(unsigned int value);

private:
	unsigned int seed;
	//	C++11 random number generators
	std::random_device rd;
	std::mt19937 gen;
};

typedef long rna_pos_type;

static std::string blankString = "";

//	Result options
static std::string ambiguousString = "__ambiguous";
static std::string noFeatureString = "__no_feature";
static std::string lowQualString = "__too_low_aQual";
static std::string notAlignedString = "__not_aligned";
static std::string notUnique = "__alignment_not_unique";

static std::string duplicateNonUnique = "__duplicateNonUnique";


//
//	The rnaPosVec class holds the set of rna fragment positions associated with one strand direction of a gene
class rnaPosVec : public std::vector<rna_pos_type>
{
public:
	//	Remove values outside the range of the annotated transcript.
	rnaPosVec & removeInvalidValues(rna_pos_type maxVal);

	//	Select a subset of the rna fragment positions
	void selectAtMost(size_t s);
};


//
//	Allows the printTsv library class to print a vector of positions
inline bool printVal(outputDataFile * f, const rnaPosVec value)
{
	return printVal(f, (const std::vector<rna_pos_type>)value);
};

//	Holds the rna-seq read position information associated with a specific gene/transcript.
class rnaSeqPositionData 
{
public:
	rnaSeqPositionData(size_t index = 0) :index(index) {};

	//	Used by the copy constructor in GeneCountData::addEntry.  Provides a more efficient way of copying the
	//	position vectors as they are just about to be discarded
	rnaSeqPositionData(size_t index, rnaPosVec && posPos, rnaPosVec && negPos) :	index(index)
	{
		swap(posPos, positions[0]);
		swap(negPos, positions[1]);
	};

	void reset() {
		positions[0].clear();
		positions[1].clear();
	}

	//	Index to the position where the associated name length and count information is held for this gene
	size_t index;

	//	The read positions themselves.  Index 0 for forward reads, 1 for reverse reads
	rnaPosVec positions[2];
};

//	This class holds and processes all of the count information associated with the set of genes
//	or transcripts
class countInfo
{
public:
	countInfo(const std::string & name,bool useForParameterEstimation) : histoGram_ind(0), 
		useForParameterEstimation(useForParameterEstimation),name(name) {};
	int histoGram_ind;
	//  Use for parameter estimation
	bool useForParameterEstimation;
	std::string name;
};

//	Hold all of the read position data in a map indexed by the gene name
class readPositionDataClass : public std::map<const stringEx, rnaSeqPositionData >
{
public:
	ADD_ITER(readGeneName, readGeneAttributes)
};




struct geneListFilenameData
{
geneListFilenameData() :start(0), finish(10000000) {};
std::string geneListFilename;
int start, finish;
};


class GeneCountData
{
	static rnaPosVec nullData;
public:

	VEC_DATA_TYPE & count(const std::string & gene);
	VEC_DATA_TYPE length(const std::string & gene);
	//	Needed if we decide the data is not name ordered and have to restart
	void reset() { for (auto & gene : readPositionData)	gene.second.reset(); };

	void addEntry(const std::string & name, bool useForParameterEstimation = false,
		VEC_DATA_TYPE length = 0, VEC_DATA_TYPE count = 0,
		rnaPosVec & posPositions = nullData, rnaPosVec & negPositions = nullData);
	void addErrorEntry(std::string name);


	std::string loadData(const std::string filename, int Nlines = -1);
	void remove_invalid_values();
	void histc (const std::vector<int> E, int maxGeneLengthForParameterEstimation);
	void transferTo(mcmcGeneData & mcmcData,size_t maxLength,int maxTotReads,int maxGeneLengthForParameterEstimation);

	bool outputGeneCounts(const std::string & filename, int detailLevel = 0, stringEx model = "");
	bool outputLandscape(const std::string & filename);
	bool outputHeatmapData(const stringEx & filename);

	//	Get read distribution 
	void getDistribution(const std::vector<int> & lengths, size_t points, std::vector<dataVec>  & counts);


	void useSelectedGenes(const geneListFilenameData & filename);

	void flatten();

	//	Names, bias, lengths and RPM data are held in a series of vectors sharing a common gene order
	std::vector<countInfo> info;
	dataVec bias;
	//	The first entry is for the raw lengths, and the second for the normalised lengths
	dataVec lengths[2];
	dataVec counts;
	dataVec RPM[2], RPKM[2], RPK[2], TPM[2];

	//	The number of genes in each range of lengths
	std::vector<double> freq;

	//	Read position data for all of the genes
	readPositionDataClass readPositionData;

	//	For data associated with reads that do not map
	readPositionDataClass errorCounts;
	std::vector<std::string> errorNames;
	dataVec errCounts;

};

#endif