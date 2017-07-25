// ***************************************************************************
// MakeFastq.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For making Fastq FIles from BamFiles
// ***************************************************************************

#ifndef MAKEFASTQ_H
#define MAKEFASTQ_H

#include <vector>
#include "api/BamReader.h"
#include "libCommon.h"
#include "stringEx.h"

using namespace std;
using namespace BamTools;

class bamRead
{
public:
	bamRead(){};
	bamRead(const BamAlignment & ba);
	void setName(const BamAlignment & ba);
	void setValues(const bamRead & br);
	void addSNP(double errorRate);
	bamRead & operator = (const BamAlignment & ba);
	void output(FILE * f);

	string readSeq,qualData;
	stringEx name; 
};

class pairedBamReads : vector<bamRead>
{
public:
	pairedBamReads() : vector<bamRead>(2)
	{
		resetNames();
	};

	void resetNames()
	{
		at(0).name = "X";
		at(1).name = "Y";
	}
	void setName(const BamAlignment & ba)
	{
		for (auto & i : This) i.setName(ba);
	}

	void addSNP(double errorRate)
	{
		for (auto & i : This) i.addSNP(errorRate);
	}
	void addRead(const BamAlignment & ba)
	{
		if (ba.IsFirstMate())
			at(0) = ba;
		else
			at(1) = ba;
	}
	bool namesMatch()
	{
		return (at(0).name == at(1).name);
	};
	void setValues(const pairedBamReads & tbr)
	{
		at(0).setValues(tbr[0]);
		at(1).setValues(tbr[1]);
	}
	void output(FILE * f1,FILE * f2)
	{
		at(0).output(f1);
		at(1).output(f2);
	}
};

class bamReadCache : public vector<pairedBamReads>
{
public:
	bamReadCache(size_t size):vector<pairedBamReads>(size){};
	void clear();
	pairedBamReads & operator[](size_t i){ return vector<pairedBamReads>::operator[](i);};
};

class MakeFastq
{
	int count;

	BamReader reader,foreignReader;
	BamAlignment ba,foreignBa;

	bool getNextAlignment();
	bool getNextForeignAlignment();
	bool getNextAlignmentCore();

	void incrementCount();

public:
	int main(int argc, char **argv);
};

#endif