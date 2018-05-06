// ***************************************************************************
// refSeqs.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 6 May 2018
// ---------------------------------------------------------------------------
// Extends featureFile class in bioinformaticsLib to support processing of 
// gff and gtf files
// ***************************************************************************

#ifndef REFSEQS_HEADER
#define REFSEQS_HEADER

#include <string>
#include <map>
#include "reference.h"
#include "hisat2Lib.h"

#include "FeatureFileEx.h"

//	When obtaining transcript sequences, this is added to either end
#define SEQ_PADDING 0

template <typename TIndexOffU> void get_sequence(
	BitPairReference<TIndexOffU>& ref,
	size_t refi,
	size_t start,
	size_t finish,
	string & sequence)
{
	sequence.resize(finish - start);
	static size_t incr = 1000;
	static uint32_t buf[(1000 + 128) / 4];
	size_t c = 0;
	ASSERT_ONLY(vector<uint32_t> destU32);
	for (size_t i = start; i < finish; i += incr) {
		size_t amt = min(incr, finish - i);
		assert_leq(amt, incr);
		int off = ref.getStretch(buf, refi, i, amt ASSERT_ONLY(, destU32));
		uint8_t *cb = ((uint8_t*)buf) + off;
		for (size_t j = 0; j < amt; j++) {
			assert_range(0, 4, (int)cb[j]);
			sequence[c++] = "ACGTN"[(int)cb[j]];
		}
	}
}


class refSeqs
{
public:
	bool get(const std::string& ebwtFileBase, std::map<std::string, geneData> & genes);

private:
	std::vector<std::string> refnames;
	std::map<std::string, size_t> nameToId;

	template <class TIndexOffU> bool getInternal(const std::string& ebwtFileBase, std::map<std::string, geneData> & genes)
	{
		vector<TIndexOffU> refLengths;
		readEbwtRefnamesAndLengths<TIndexOffU>(ebwtFileBase, refnames, refLengths);

		for (size_t i = 0; i < refnames.size(); i++)
		{
			string id;
			parser(refnames[i], " ", id);
			nameToId.emplace(id, i);
		}


		BitPairReference<TIndexOffU> ref(
			ebwtFileBase, // input basename
			false,                // sanity-check reference
			false,              // be talkative
			false);             // be talkative at startup
		assert_eq(ref.numRefs(), refnames.size());

		size_t geneCount = 0;

		progMessage("Getting sequences from reads");

		for (auto & gene : genes)
		{
			if ((++geneCount % 2000) == 0)
				progMessage(geneCount, " genes processed");
			auto i = nameToId.find(gene.second.chromosome);
			featureRegionList & regions = gene.second.regions;
			gene.second.length = 0;
			if (i != nameToId.end())
			{
				size_t refId = i->second;
				long long int start = regions[0]->start - 2;
				if (start < 1)
				{
					gene.second.priorSeq.assign(SEQ_PADDING, 'N');
				}
				else if (start < SEQ_PADDING)
				{
					get_sequence<TIndexOffU>(ref, refId, 1, start, gene.second.priorSeq);
					gene.second.priorSeq = string(SEQ_PADDING - gene.second.priorSeq.size(), 'N') + gene.second.priorSeq;
				}
				else
					get_sequence<TIndexOffU>(ref, refId, start - SEQ_PADDING, start, gene.second.priorSeq);

				size_t end = regions.back()->finish + 1;
				if (end > refLengths[refId])
				{
					gene.second.postSeq.assign(SEQ_PADDING, 'N');
				}
				else if (end > (refLengths[refId] - SEQ_PADDING))
				{
					get_sequence<TIndexOffU>(ref, refId, end, refLengths[refId], gene.second.postSeq);
					gene.second.postSeq += string(SEQ_PADDING - gene.second.postSeq.size(), 'N');
				}
				else
					get_sequence<TIndexOffU>(ref, refId, end, end + SEQ_PADDING, gene.second.postSeq);

				for (int j = 0; j < regions.size(); j++)
				{
					gene.second.mRNAtoSeq.emplace(regions[j]->RNAstart, regions[j]->start - 1);
					get_sequence<TIndexOffU>(ref, refId, regions[j]->start - 1, regions[j]->finish, regions[j]->sequence);
					gene.second.length += (regions[j]->finish - regions[j]->start + 1);
				}
			}
			else
			{
				gene.second.priorSeq.assign(SEQ_PADDING, 'N');
				gene.second.postSeq.assign(SEQ_PADDING, 'N');

				for (int j = 0; j < regions.size(); j++)
				{
					gene.second.mRNAtoSeq.emplace(regions[j]->RNAstart, regions[j]->start - 1);
					regions[j]->sequence.assign(regions[j]->finish - regions[j]->start + 1, 'N');
					gene.second.length += (regions[j]->finish - regions[j]->start + 1);
				}
			}
		}
		return true;
	}
};

#endif