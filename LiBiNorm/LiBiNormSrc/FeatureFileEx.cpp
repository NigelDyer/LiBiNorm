// ***************************************************************************
// FeatureFileEx.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Extends featureFile class in bioinformaticsLib to support processing of 
// gff and gtf files
// ***************************************************************************#include "Options.h"
#include "FeatureFileEx.h"

using namespace std;

// Check to see if the feature overlaps a region as defined by a segment.  If it does then add the
//	information to the overlap list.  A strict overlap is when the feature fully overlaps the segment
void featureRegion::checkOverlap(const region & segment,vector<featureOverlap> & overlapList) const
{
	rna_pos_type RNAstartPos = max<long>(segment.start - start + RNAstart,0);
	rna_pos_type RNAendPos = segment.end-start + RNAstart;

	//	strict fit
	if ((segment.start >= start) && (segment.end <= finish))
	{
		overlapList.emplace_back(segment.start,segment.end, RNAstartPos, RNAendPos,true,name,type);
	}
	//	non-strict fits
	else if ((segment.start <= finish) && (segment.start >= start)) 
	{
		overlapList.emplace_back(segment.start,finish, RNAstartPos, RNAendPos, false,name,type);
	}
	else if ((segment.end <= finish) && (segment.end >= start))
	{
		overlapList.emplace_back(start,segment.end, RNAstartPos, RNAendPos, false,name,type);
	}
	else if ((segment.start <= start) && (segment.end >= finish))
	{
		overlapList.emplace_back(start,finish, RNAstartPos, RNAendPos,false,name,type);
	}
}

//	Two constructors.   The first is an R-value copy constructor for efficiency, and the second constructs a feature region
//	from the associated values
featureRegion::featureRegion(featureRegion && gtf) : start(gtf.start), finish(gtf.finish), RNAstart(gtf.RNAstart), name(std::move(gtf.name)), type(std::move(gtf.type)),
bioType(std::move(gtf.bioType)), strand(gtf.strand), overlaps(gtf.overlaps)
{};

featureRegion::featureRegion(size_t start, size_t finish, const std::string & name, char strand, const std::string & type, const std::string & bioType) :
	start(start), finish(finish),  RNAstart(0), name(name), type(type), bioType(bioType),strand(strand)
{};


void featureFileEx::index(GeneCountData & geneCounts, bool useStrand)
{
	bool usingPreselectedGenes = false;
	//	If there are already entries in the geneCounts data at this stage it is because
	//	we have preloaded them with a set of genes/transcripts that we are specifically
	//	interested in.  At this point we then get rid of the rest
	if (geneCounts.readPositionData.size() > 1)
	{
		usingPreselectedGenes = true;
		for (entryMapClass::Pair chrom : entryMap)
		{
			//	First get rid of entries associated with genes that we are not interested in 
			for (featureFile::featureMap::Iterator i = chrom.data().begin(); i != chrom.data().end();)
			{

				string prefix(i.feature().name.substr(0, i.feature().name.find('.')));
				readPositionDataClass::Iterator j = geneCounts.readPositionData.lower_bound(prefix);
				//	Need to do this to cater for genes/transcripts that have changed release/version 
				if (!j.readGeneName().startsWith(prefix))
				{
					auto k = i++;
					chrom.data().erase(k);
				}
				else
				{
					if (i.feature().name != j.readGeneName())
					{
						//	We now have a different version of the transcript/gene to that in the reference list, so rename
						//	our list to match the new reference genome
						geneCounts.readPositionData.emplace(i.feature().name, j.readGeneAttributes().index);
						geneCounts.info[j.readGeneAttributes().index].name = i.feature().name;
						geneCounts.readPositionData.erase(j);
					}
					i++;
				}

			}

		}
	}

	for (entryMapClass::Pair chrom : entryMap)
	{
		//	In each chromosome go through all of the regions to see what regions can be amalgamated
		featureRegion::chromosomeFeatureData & thisChromData = genomeGtfData[chrom.name()];
		for (featureFile::featureMap::Iterator i = chrom.data().begin(); i != chrom.data().end();i++)
		{

			size_t finish = i.feature().finish;
			setEx<string> type{{i.feature().type}};
			for (featureFile::featureMap::Iterator j = next(i,1);(j != chrom.data().end()) && (j.position() <= (finish + 1));)
			{
				featureFile::featureMap::Iterator k = j++;
				if (k != chrom.data().end())
				{
					//	It turns out that the Yam1 gene is defined on both strands, so need to check strands
					if ((i.feature().type == k.feature().type) && 
						(i.feature().name == k.feature().name) && 
						(i.feature().strand == k.feature().strand))
					{
						//	If the second region extends beyond the firat then make the region longer and add
						//	the second type to the list of types associated with the region
						if (k.feature().finish > finish)
						{
							finish = k.feature().finish;
							type.add(k.feature().type);
							chrom.data().erase(k);
						}
						//	If the second region is the same length or shorter then just add the type associated with the new region
						else if (k.feature().finish <= finish)
						{
							type.add(k.feature().type);
							chrom.data().erase(k);
						}
					}
				}
			}
			static const string nullString;
			thisChromData.emplace(i.position(),featureRegion(i.feature().start,finish,i.feature().name,
				i.feature().strand,i.feature().type , i.feature().biotype));
		}

		//	And now for each region find the regions that it overlaps and produce overlap list:  The list of all regions that start before this region has ended.
		//	Then select the first of the regions, which will be used as the starting point for checking for region overlaps
		chromosomeEndIndexMap & thisChromEndMap = genomeEndIndex[chrom.name()];
		
		//	We are using a temporary map of the address of the features in this chromosome.
		//	We can use address only because the entries will not be moved during this process
		//	For each of the features we create a map of iterators pointing to other overlapping features, 
		//	the mmap being indexed by the position of the start of the region
		class overlapMapType : public map<const void *, chromosomePositionIndexMap >
		{
		public:
			ADD_ITER(feature,positionData)
		} overlapMap;

		for (featureRegion::chromosomeFeatureData::Iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			//	Take the opportunity to produce a map of all the genes for holding counts
			//	The default is that it will be used for parameter estimation
			geneCounts.addEntry(i.feature().name,true);

			//	And a parallel map of the ends of all of the featureRegions/
			thisChromEndMap.emplace(i.feature().finish,i);

			//	Starting from the next region, find all of the subsequent regions which start before this region ends.  
			//	In each case add the subsequent region to the list of overlaps.  And also add this region to the list of 
			//	overlaps
			for (featureRegion::chromosomeFeatureData::Iterator j = next(i, 1); (j != thisChromData.end()) && (j.position() < i.feature().finish); j++)
			{
				overlapMap[&j.feature()].emplace(i.position(), i);
//				overlapMap[&i->second].emplace(j->first, j);
			}
		}
		for (featureRegion::chromosomeFeatureData::Iterator i = thisChromData.begin(); i != thisChromData.end();i++)
		{
			size_t index = geneCounts.readPositionData.at(i.feature().name).index;

			//	Use this for accumulating length information
			VEC_DATA_TYPE & length = geneCounts.lengths[0].at(index);
			overlapMapType::Iterator j = overlapMap.find(&i.feature());

			if (j == overlapMap.end())
			{
				//	The overlaps pointer points to self, indicating that there is no overlap
				i.feature().overlaps = i;
				//	Add featureRegion to the list of regions associated with the gene
				genes[i.feature().name].addRegion(&i.feature(), false,length,chrom.name());

			}
			else
			{
				i.feature().overlaps = j.positionData().Begin().featureDataIterator();
				genes[i.feature().name].addRegion(&i.feature(), true,length,chrom.name());
				//	If the genes overlap then we dont use them if either
				//	   a) the reads are unstranded
				//	   b) the genes are on the same strand
				if (!usingPreselectedGenes && 
					(!useStrand || (i.feature().strand == j.positionData().Begin().featureDataIterator().feature().strand)))
				{
					geneCounts.info[index].useForParameterEstimation = false;
					size_t index2 = geneCounts.readPositionData.at(j.positionData().Begin().featureDataIterator().feature().name).index;
					geneCounts.info[index2].useForParameterEstimation = false;
				}
			}
		}
	}
	geneCounts.addErrorEntry(noFeatureString);
	geneCounts.addErrorEntry(ambiguousString);
	geneCounts.addErrorEntry(lowQualString);
	geneCounts.addErrorEntry(notAlignedString);
	geneCounts.addErrorEntry(notUnique);
}

//	Prints a list of all the features, ordered by chromosome and then 

bool featureFileEx::outputChromData(const string & filename, const GeneCountData & geneCounts)
{
	TsvFile output;
	output.open(filename);

	if (!output.is_open())
	{
		progMessage("Unable to open ", filename," for outputting genome data");
		return false;
	}

	for(auto & i : genomeGtfData)
	{
		for (featureRegion::chromosomeFeatureData::Pair j : i.second)
		{
			size_t index = geneCounts.readPositionData.at(j.feature.name).index;
			size_t countF = geneCounts.readPositionData.at(j.feature.name).positions[0].size();
			size_t countR = geneCounts.readPositionData.at(j.feature.name).positions[1].size();
			bool beingUsed = geneCounts.info[index].useForParameterEstimation;
			int length = geneCounts.lengths[0][index];
			output.printEnd(i.first, _s("chr",i.first, ":", j.feature.start, "-", j.feature.finish), j.feature.RNAstart, j.feature.strand,
				j.feature.name, j.feature.type, j.feature.bioType,beingUsed, length,countF, countR);
		}
	}
	return true;
}

