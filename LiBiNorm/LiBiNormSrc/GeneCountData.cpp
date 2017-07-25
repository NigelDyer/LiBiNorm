// ***************************************************************************
// GeneCountData.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For processing count information associated with genes
// ***************************************************************************#include <fstream>

#include "containerEx.h"
#include "Options.h"
#include "GeneCountData.h"

using namespace std;

//  Remove reads that, apparently, start before the beginning or after the end of the gene.
//	The first version reproduces the less efficient algorithm that was used in th eoriginal matlab
//	code
rnaPosVec & rnaPosVec::removeInvalidValues(rna_pos_type maxVal)
{
#ifdef REPRODUCE_MATLAB
	// Reproduce the matlab code so that results can be compared
	size_t offset = 0;
	for (size_t i = 0; i + offset < size();)
	{
		if ((at(i + offset) < 0) || (at(i + offset) >= maxVal))
			offset++;
		else
		{
			if (offset)
				at(i) = at(i + offset);
			i++;
		}
	}
	if (offset)
		resize(size() - offset);
#else
	//	Sort in place for maximum efficiency, if we find an invalid value, replace with one from the end;
	//	Note that if we swap with a value from the end we have to check it as well to see if it is invalid
	iterator i = begin(), j = end();
	while (i != j)
	{
		if ((*i < 0) || (*i >= maxVal))
			std::swap(*i, *--j);
		else
			i++;
	}
	resize(j - begin());
#endif
	return *this;
}


/*
//  Returns a random number, uniformly distributed between 0 and a.
int intRand(int max)
{
	// Static members generate distribution in range 0 to 1.   This is scaled by input parameter for each call
	static std::random_device rd;
#ifdef SELECT_READS_SEED
	static std::mt19937 gen(SELECT_READS_SEED);
#else
	static std::mt19937 gen(rd());
#endif

	static std::uniform_int_distribution<> dis(0, 65535);

	return dis(gen) % max;
}

*/

intRandClass intRand;

//	Test mode as set by REPRODUCE_MATLAB takes the first N samples rather than randomly picks samples,
//	and uses the same algorithm as the MATLAB code for excluding invalid calues.   Used for comparing the two outputs
void rnaPosVec::selectAtMost(size_t s)
{
	if (s > size())
		return;
	//For reproducibility in test mode just use the first N values
#ifndef REPRODUCE_MATLAB
	//	Swap the first s entries with the entry at some other position, then resize to just have the s entries
	iterator i = begin();
	for (size_t j = 0; j < s; j++)
	{
		std::swap(*(i++), *(begin() + intRand.value(size())));
	}
#endif
	resize(s);
}


rnaPosVec GeneCountData::nullData;


//	Returns the address of the counter associated with a gene or error condition
//	such that the value can be incremened when an associated error if found

VEC_DATA_TYPE & GeneCountData::count(const std::string & gene)
{
	static VEC_DATA_TYPE dummy;
	readPositionDataClass::Iterator i = readPositionData.Find(gene);
	if (i == readPositionData.end())
	{
		i = errorCounts.Find(gene);
		if (i == errorCounts.end())
		{
			_ASSERT_EXPR(false, "Looking for gene that was not in the reference genome");
			return dummy;
		}
		return errCounts.at(i.readGeneAttributes().index);
	}
	return counts.at(i.readGeneAttributes().index);
}

//	Returns the actual length of a gene (rather than the normalised length) given the name  
VEC_DATA_TYPE GeneCountData::length(const std::string & gene)
{
	readPositionDataClass::Iterator i = readPositionData.find(gene);
	if (i == readPositionData.end())
	{
		return 0;
	}
	return lengths[0].at(i.readGeneAttributes().index);
}

//	Reads in a file contain a list of gene names.  Only these genes will then be used
//	in the subsequent analysis
void GeneCountData::useSelectedGenes(const geneListFilenameData & filenameData)
{
	ifstream file;
	file.open(filenameData.geneListFilename);

	if (!file.is_open())
	{
		progMessage("Unable to read gene list from ", filenameData.geneListFilename);
		return;
	}

	stringEx line, gene;
	int i = 0;
	while (!file.eof())
	{
		std::getline(file, line);
		parser(line, " \n\r",gene);
		if ((i >= filenameData.start) && (i <= filenameData.finish) && (gene))
		{
			addEntry(gene,true);
		}
		i++;
	};
}

//	Print the results to varying levels of detail
//	0 = htseq-count compatible
//	1 = basic normalised results
//	2 = As 1 but with header
//	3 = As 2 but more columns

bool GeneCountData::outputGeneCounts(const string & filename, int detailLevel, stringEx model)
{
	TsvFile output;

	model.replace(" ", "");

	if (!output.open(filename))
		return false;

	if ((bias.size()) && (detailLevel > 0))
	{
		lengths[1] = lengths[0] * bias;

		// Definitions taken from http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
		for (size_t i = 0; i < ((bias.size()) ? 2 : 1); i++)
		{
			VEC_DATA_TYPE scalingFactor = sum(counts) / 1000000;

			RPM[i] = counts / scalingFactor;
			RPKM[i] = RPM[i] / lengths[i];

			RPK[i] = counts / lengths[i];
			scalingFactor = sum(RPK[i]) / 1000000;
			TPM[i] = RPK[i] / scalingFactor;
		}

		switch (detailLevel)
		{
		case 1:
			output.print("Gene", "count", "length", "Bias", _s("TPM_", model));
			break;
		case 2:
			output.print("", "", "RNA", "Bias:", model + " +");
			output.print("Gene", "count", "length", model, "TPM");
			output.print();
			break;
		case 3:
			output.print("", "", "RNA", "Bias:", model + " +", model + " +", model + " +", model + " +", "Raw", "Raw", "Raw", "Raw");
			output.print("Gene", "count", "length", model, "RPM", "RPKM", "RPK", "TPM", "RPM", "RPKM", "RPK", "TPM");
			output.print();
			break;
		}

		//	Dont start at 0 as 0 is the reference for normalisation
		for (size_t i = 1; i < info.size(); i++)
		{
			output.printStart(info[i].name, counts[i]);

			switch (detailLevel)
			{
			case 1:
			case 2:
				output.printMiddle(lengths[0][i], bias[i], TPM[1][i]);
				break;
			case 3:
				output.printMiddle(lengths[0][i], bias[i],
					RPM[1][i], RPKM[1][i], RPK[1][i], TPM[1][i],
					RPM[0][i], RPKM[0][i], RPK[0][i], TPM[0][i]);
				break;
			}
			output.printEnd();
		}
	}
	else
	{
		//	We are just printing the counts in htseq-count mode, which means that they
		//	need to be in name order, and not in the order they are found 
		for (auto i : readPositionData)
		{
			if (i.second.index != 0)
				output.print(i.first, (int)counts[i.second.index]);
		}
	}

	//	Print the error counts at the end
	if (detailLevel != 1)
	{
		for (size_t i = 0; i < errorNames.size(); i++)
			output.print(errorNames[i], (int)errCounts[i]);
	}

	return true;
}


bool GeneCountData::outputLandscape(const string & filename)
{
	TsvFile output;

	if (!output.open(filename))
		exitFail("Unable to open ", filename, " for landscape data");
#ifdef LANDSCAPE_FORMAT_1
	//	This is the original format
	for (size_t i = 1;i < info.size();i++)
	{
		long count = counts[i];
		if (/*(count > 0) && */(info[i].useForParameterEstimation))
		{
			long len = lengths[0][i];
			string & name = info[i].name;
			string c;
#ifdef COUNT_IN_LANDSCAPE
			c = _s(":", count);
#endif
			output.print(name, _s(len, c, " plus"), readPositionData[name].positions[0]);
			output.print(name, _s(len, c, " minus"), readPositionData[name].positions[1]);
		}
	}
#else
	//	This format includes an additional header line, and an extra line per gene/transcript
	//	that contains the meta information such as length, number of reads (redundant but
	//	useful when viewing the file and whether it is 
	output.print("Landscape file","Format",2);
	for (size_t i = 1; i < info.size(); i++)
	{
		long count = counts[i];
		long len = lengths[0][i];
		string & name = info[i].name;
		output.print(name, len, count, info[i].useForParameterEstimation?"Y":"N");
		output.print(name, "plus", readPositionData[name].positions[0]);
		output.print(name, "minus", readPositionData[name].positions[1]);
	}

#endif

	return true;
}

//	Enable this to produce output based on peak rather than average values
// #define PEAKVALUES

//	Produces a datafile that can be used by LiBiNormPlot.R to produce a heat map of the read distribution
bool GeneCountData::outputHeatmapData(const stringEx & filename)
{
#ifdef MATHMATICA_FILE
	CsvFile mathematicaFile;
	string mFilename = filename.replaceSuffix(".mm.txt");
	if (!mathematicaFile.open(mFilename))
		exitFail("Unable to open ", mFilename, " for bias data");
	fprintf(mathematicaFile.fout, "{");
	bool first = true;
#endif
	TsvFile tsvFile;
	if (!tsvFile.open(filename))
	{
		progMessage("Unable to open ", filename, " for bias data");
		return false;
	}
	multimap<VEC_DATA_TYPE, size_t> orderedLengths;
	//	Create a map of all the genes, ordered by length.  Ignore genes with no reads
	for (size_t i = 1; i < info.size(); i++)
	{
		auto & data = readPositionData[info[i].name];
		size_t Nreads = data.positions[0].size() + data.positions[1].size();
		if (Nreads > READ_THRESHOLD)
			orderedLengths.emplace(lengths[0][i], i);
	}

	dataArray allBins(orderedLengths.size(), N_BIAS_BINS);
	size_t dataCount = 0;
	for (auto i = orderedLengths.begin();i != orderedLengths.end();i++)
	{
		//	Get info for the gene
		dataVec & counts = allBins[dataCount++];
		string & name = info[i->second].name;
		auto & data = readPositionData[name];
		VEC_DATA_TYPE length = i->first;

		//And then setup the bins
		dataVec E(N_BIAS_BINS);
		for (int j = 0; j < N_BIAS_BINS; j++)
			E[j] = (length * j)/ N_BIAS_BINS;

		size_t Nreads = data.positions[0].size() + data.positions[1].size();
		//	Put all the reads into bins
		for (size_t strand = 0; strand < 2; strand++)
		{
			rnaPosVec & readPositions = data.positions[strand];
			for (auto v : readPositions)
			{
				size_t l = 0;
				size_t h = E.size() - 1;
				size_t k = l;
				while ((h - l) > 1)
				{
					k = (h + l) / 2;
					if (v < E[k])
						h = k;
					else
						l = k;
				}
				if (v >= E[h])
					k = h;
				else
					k = l;
				counts[k]++;
			}
		}
		//	and normalise.  The mathematica code uses counts normalised by reads
		counts /= Nreads;
#ifdef MATHMATICA_FILE
		if (first)
			first = false;
		else
			fprintf(mathematicaFile.fout, ",");

		fprintf(mathematicaFile.fout, "{");
		mathematicaFile.printStart(fmt("%f", counts));
		fprintf(mathematicaFile.fout, "}");
#endif

/*		//	And then scale the bins so that they go between 0 and 1
		VEC_DATA_TYPE min = 10, max = 0;
		for (auto n : counts)
		{
			if (n < min) min = n;
			if (n > max) max = n;
		}
		for (auto & n : counts)
		{
			n = (n - min) / (max - min);
		}*/
	}
#ifdef MATHMATICA_FILE
	fprintf(mathematicaFile.fout, "}");
#endif


#ifdef N_BIAS_GENE_SEGMENTS
	//	And now consolidate the data down to N_BIAS_GENE_SEGMENTS (500) rows by combining genes
	dataArray consolidatedBins(N_BIAS_GENE_SEGMENTS, N_BIAS_BINS);

	size_t lastPos = 0;
	size_t section_count = 0;
	for (size_t i = 0; i < allBins.size(); i++)
	{
		size_t currPos = (i * N_BIAS_GENE_SEGMENTS)/ allBins.size();
		if (currPos != lastPos)
		{
#ifndef PEAKVALUES
			//	Normalise by the number of genes in this section
			consolidatedBins[lastPos] /= section_count;

			//	This does infill if there are fewer than N_BIAS_GENE_SEGMENTS genes
			for (size_t i = lastPos + 1; i < currPos; i++)
				consolidatedBins[i] = consolidatedBins[lastPos];
#endif
			section_count = 0;
			lastPos = currPos;
		}
#ifdef PEAKVALUES
		for (size_t j = 0; j < bins; j++)
			if (allBins[i][j] > consolidatedBins[currPos][j])
				consolidatedBins[currPos][j] = allBins[i][j];
#else
		consolidatedBins[currPos] += allBins[i];
#endif
		section_count++;
	}
#ifndef PEAKVALUES
	consolidatedBins[lastPos] /= section_count;

	//	Clip values above MAX_VALUE
	for (auto & i : consolidatedBins)
		i.clipMax(MAX_VALUE);

#endif
	//	Output in reverse order so that it is plotted correctly by R
	for ( int i = consolidatedBins.size()-1; i >= 0;i--)
		tsvFile.print(fmt("%f", consolidatedBins[i]));
#else
	//	This is the code for if we are not consolidating but outputting all of the genes
	for (int i = allBins.size() - 1; i >= 0; i--)
		tsvFile.print(fmt("%f", allBins[i]));
#endif

	return true;
}

//	Find the number of read position values in 'this' (a vector) that sits within each bin of the
//	histogram defined by E.   The results go into the 'freq'vector
//	This is used as part of the liklyhood calculations
//	Only include genes/transcripts that are OK for use during the parameter estimation phase,
//	ie do not overlap other genes, and are less than the current length threshold.
void GeneCountData::histc(const vector<int> E, int maxGeneLengthForParameterEstimation)
{
	freq.assign(E.size(), 0);

	//	Increment from 1 because the first entry is the reference length
	for (size_t i = 1; i < info.size(); i++)
	{
		if ((info[i].useForParameterEstimation) && (lengths[0][i] <= maxGeneLengthForParameterEstimation))
		{
			double v = lengths[0][i];
			size_t l = 0;
			size_t h = E.size() - 1;
			size_t k = l;
			while ((h - l) > 1)
			{
				k = (h + l) / 2;
				if (v < E[k])
					h = k;
				else
					l = k;
			}
			if (v == E[h])
				k = h;
			else
				k = l;
			info[i].histoGram_ind = k;
			freq[k]++;
		}
	}
}

void GeneCountData::remove_invalid_values()
{
	for (size_t i = 0; i < info.size(); i++)
	{
		for (size_t j = 0; j < 2; j++)
			readPositionData[info[i].name].positions[j].removeInvalidValues(lengths[0][i]);
	}
}


void GeneCountData::addEntry(const string & name, bool useForParameterEstimation,
	VEC_DATA_TYPE length, VEC_DATA_TYPE count,
	rnaPosVec & posPositions, rnaPosVec & negPositions)
{
	auto geneData = readPositionData.find(name);
	if (geneData == readPositionData.end())
	{
		counts.push_back(count);
		info.push_back(countInfo(name, useForParameterEstimation));
		lengths[0].push_back(length);
		readPositionData.emplace(name, geneReadData(counts.size() - 1, move(posPositions), move(negPositions)));
	}
}

void GeneCountData::addErrorEntry(string name)
{
	auto geneData = errorCounts.find(name);
	if (geneData == errorCounts.end())
	{
		errCounts.push_back(0);
		errorNames.push_back(name);
		errorCounts.emplace(name, errCounts.size() - 1);
	}
}


string GeneCountData::loadData(const string filename, int Ngenes)
{
	std::ifstream f;
	int format = -1;
	f.open(filename);
	if (!f.is_open())
		exitFail("Unable to open file ", filename);

	//The reference value
	addEntry("reference",false, DEFAULT_NORMALISATION_GENE_LENGTH);


	stringEx buffer;
	string gene, direction;
	int a = 1;
	stringEx lastGene;
	while (!f.eof() & (Ngenes-- != 0))
	{
		//	Need to check for gene and length consistency and that the gene has not appeared before

		getline(f, buffer);
		if (format == -1)
		{
			string dummy;
			if (buffer.startsWith("Landscape file"))
			{
				parseTsv(buffer, dummy, dummy, format);
				if ((format < 2) || (format > 2))
					exitFail("Unknown landscape file format");
				getline(f, buffer);
			}
			else
				format = 1;
		}
		string gene,gene1,gene2, n1, n2;
		if (buffer.size())
		{
			rnaPosVec posPositions, negPositions;
			int length, count;
			string use("Y");

			switch (format)
			{
			case 1:
				//	Need to case the posPositions to the undelying type in order to pick up the
				//	librray code for parsing a vector of things
				parseTsv(buffer, gene, n1, " ", direction, (std::vector<rna_pos_type> &)posPositions);
				getline(f, buffer);
				parseTsv(buffer, gene1, n2, " ", direction, (std::vector<rna_pos_type> &)negPositions);
				if ((gene1 != gene) || (direction != "minus"))
					exitFail("Landscape file format error: ", buffer);
				if (strchr(n1.c_str(), ':') == NULL)
				{
					length = atoi(n1.c_str());
					count = posPositions.size() + negPositions.size();
				}
				else
					parser(n1, ":", length, count);
				break;
			case 2:
				parseTsv(buffer, gene, length,count,use);
				getline(f, buffer);
				parseTsv(buffer, gene1, direction, (std::vector<rna_pos_type> &)posPositions);
				if ((gene1 != gene) || (direction != "plus")) 
					exitFail("Landscape file format error: ", buffer);
				getline(f, buffer);
				parseTsv(buffer, gene2, direction, (std::vector<rna_pos_type> &)negPositions);
				if ((gene2 != gene) || (direction != "minus"))
					exitFail("Landscape file format error: ", buffer);
				break;
			}


				addEntry(gene, use == "Y",length, count, posPositions, negPositions);
			a++;
		}
		lastGene = gene;
	}

	if(lastGene)
		progMessage("Last gene = ", lastGene);
	return lastGene;
}

//	Transfers information for up to maxLength reads from up to Ngenes genes or transcripts
//	into the form which can be used by the mcmc chain
void GeneCountData::transferTo(mcmcGeneData & mcmcData, size_t maxLength, int maxTotReads, int maxGeneLengthForParameterEstimation)
{
	vectorEx<int> bins{ { 0,300 } };
	for (size_t i = 500; i <= 10000; i += 500)
		bins.push_back(i);
	bins.add(11000, 12000, 15000, 30000);

	histc(bins, maxGeneLengthForParameterEstimation);

	freq[0] = freq[0] * 2;
	freq[1] = freq[1] * 2;
	freq[21] = freq[21] / 2;
	freq[22] = freq[22] / 2;
	freq[23] = freq[23] / 6;
	freq[24] = freq[24] / 6;

	size_t geneIndex = 0;
	mcmcData.geneLengths.resize(readPositionData.size());
	mcmcData.geneFrequencies.resize(readPositionData.size());

	srand((unsigned)time(NULL));

	int Nreads = 0;

	//	Start at one because we dont include the reference gene (except there are no reads
	//	in the reference gene so this makes no difference
	for (size_t i = 1; i < info.size(); i++)
	{
		//	Only include genes/transcripts that are OK for use during the parameter estimation phase,
		//	ie do not overlap other genes, and are less than the current length threshold.
		if ((info[i].useForParameterEstimation) && (lengths[0][i] <= maxGeneLengthForParameterEstimation))
		{
			//	For the forward and the reverse counts
			for (size_t j = 0; j < 2; j++)
			{
				rnaPosVec & positions = readPositionData.at(info[i].name).positions[j];
				positions.selectAtMost(maxLength);

				//	fragData contains the count 
				size_t len = min(maxLength, positions.size());
				mcmcData.fragData.append(positions,len);

				mcmcData.geneIndex.insert(mcmcData.geneIndex.end(),len, geneIndex);
				Nreads += positions.size();
			}
			mcmcData.geneLengths[geneIndex] = lengths[0][i];
			mcmcData.geneFrequencies[geneIndex] = freq[info[i].histoGram_ind];

			geneIndex++;
			if ((maxTotReads) && (Nreads > maxTotReads))
				break;
		}
	}
	mcmcData.geneLengths.resize(geneIndex);
	mcmcData.geneFrequencies.resize(geneIndex);

	optMessage(Nreads, " reads used for parameter determination");
}


