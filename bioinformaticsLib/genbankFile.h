// ***************************************************************************
// genbankFile.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Code for processing genbank files
// ***************************************************************************

#ifndef GENBANKFILE_H
#define GENBANKFILE_H

#include <fstream>
#include "containerEx.h"
#include "stringEx.h"
#include "printEx.h"

//	The values for the CDS start and finish in the genbank files are either integers or qualified integers, ie
//	123 or <123 or >123.    This class stores the number and the qualifier.  PrintVal is used by PrintTsv
struct qualInt
{
	qualInt(): value(0),qualifier(' '){};
	std::string asString();
	int operator - (const qualInt & a){return value-a.value;};
	const qualInt & operator = (int a){value=a;qualifier = ' ';return *this;}; 
	operator int() {return value;};
	int value;
	char qualifier;
};

class  extendableQuotedString  : public std::string
{
public:
	extendableQuotedString & operator = (const std::string & S) { std::string::operator = (S); return *this;};
};

class extendableRemainderOfLine : public stringEx 
{
public:
	extendableRemainderOfLine(const std::string &  a = "") :stringEx(a)
	{
	}
	extendableRemainderOfLine & operator = (const std::string &  a)
	{
		std::string::operator=(a);
		return *this;
	}
};

namespace parserInternal
{
	void parseval(const char *& start,qualInt & value,size_t & len);
//	Finds the next string, which may be in quotes and adds it to the string in size
	void parseval(const char *& start,extendableQuotedString & value,size_t & len);
//	Adds the text in the remainder of the line to value
	void parseval(const char * start,extendableRemainderOfLine & value,size_t & len);
}

int getVariant(const std::string & description);

inline bool printVal(outputDataFile * f,const qualInt & value)
{
	if (value.qualifier != ' ')
		fprintf(f->fout,"%c%i",value.qualifier,value.value);
	else
		fprintf(f->fout,"%i",value.value);
	return true;
};


class genbankFile;

class featureSection : public std::multimap<std::string,std::string>
{
public:
	featureSection():linkedCDS(-1),variant(0){};
	featureSection(genbankFile * f,const std::string & type,const char *& c);

	const std::string & identifier()
	{
		const static std::string unnamed("unnamed"); 
		const_iterator fi = find("transcript_id");
		if (fi != end())
			return (fi->second);
		fi = find("gene");
		if (fi != end())
			return (fi->second);
		return unnamed;
	};

	void addSuffix(const std::string & S)
	{
		iterator fi = find("transcript_id");
		if (fi != end())
			fi->second += S;
		fi = find("gene");
		if (fi != end())
			fi->second += S;
	}

	const std::string & operator [] (const std::string & index) const
	{
		const static std::string nullString; 
		const_iterator fmi = find(index);
		return ((fmi == end())?nullString:fmi->second);
	};
	std::string type;
	qualInt start,finish;
	int linkedCDS;
	int variant;
};

class featureMap: public std::vector<featureSection>
{
}; // CDS, gene,source,ncRNA;


class orfInfo
{
public:
	orfInfo():start(0),finish(0),finish2(0),score(0){};
	int start,finish,finish2;
	double score;
};

//	And then a support code for printing orf data to a tab separated variable file
bool printVal(outputDataFile * f,const orfInfo & value);


class nucleotideSequence : public std::string
{
public:
	void getSecondORFs(size_t orf1Start,int geneEnd,std::vector<orfInfo> & orfs);
	size_t getSecondStopCodon(size_t orf1Start);
	bool seqMatch(size_t offset,const char * seq);
};

class geneInfo
{
	const std::string & getVal(const char * str1,const char * str2);
public:
	featureMap features;
	std::string type;

	const char * readNext(genbankFile * f);

	const char * geneName();
	const std::string & ncRNA_class();
};

class genbankFile;

class genbankLocus 
{
	friend class featureSection;
public:
	bool readNext(genbankFile * f);
	bool readOriginSeq(genbankFile * f);
	const std::string & mol_type();
	const char * version()
	{
		if (VERSION == "")
			return LOCUS_id.c_str();
		return VERSION.c_str();
	}
	featureSection source; // CDS, gene,source,ncRNA;
	stringEx LOCUS_id, VERSION, GI;
	extendableRemainderOfLine SOURCE,DEFINITION;
	int LOCUS_length,variant; 
	nucleotideSequence nucSeq;
};


class genbankFile
{
	std::ifstream f;
	std::string buffer;
public:
	const char * line()
	{
		return buffer.c_str();
	}
	const char * getLine()
	{
		if (!f.eof())
		{
			if (getline(f,buffer))
				return buffer.c_str();
		}
		buffer.clear();
		return nullptr;
	};



	bool open(const char * filename);
	bool open(const std::string & filename){ return open(filename.c_str());};
	void rewind() {f.clear();f.seekg(0,std::ios::beg);};

	bool getNextLocus(genbankLocus & locus);
	bool getNextGene(geneInfo & gene);
};

/**************************************************************
	For analysing a set of reads associated with the genome

*/

class genbankInfoParser 
{
	friend struct genbankInfo;
	enum columnTypes {
	null,geneIdPos,geneNamePos,chromosomePos,startPos,endPos,lengthPos,typePos};
	std::vector<columnTypes> columns;
public:
	genbankInfoParser(const stringEx & header);
};

struct genbankInfo
{
	genbankInfo():length(0),orf2offset(0),orf2start(0),orf2finish(0),orf2finish2(0),orf2length(0),reads(0),validReads(0),orf2score(0)
	{
	};

	genbankInfo(const char * initialiserString);
	genbankInfo(const std::string &  initialiserString,const genbankInfoParser & info,std::string & id);
	int length,orf2offset,orf2start,orf2finish,orf2finish2,orf2length,reads,validReads;
	double orf2score;

	std::string type,gene,ncRNAtype,definition;
	std::vector<stringEx> variants;
	std::string chromosome;
	qualInt geneStart,CDSstart,CDSfinish;
	int finish2,variant;
};

class genInfoFile : public std::map<stringEx,genbankInfo>
{
public:
	class chromosomePositions : public map<int,std::string> {};
	map<std::string,chromosomePositions> chromosomeInfo; 
	bool load(const std::string & filename);
};

/**************************************************************
	For processing sets of genBank transcript data as retreived from a genbank file

*/

struct genbankTranscriptInfo
{
	genbankTranscriptInfo():length(0),reads(0),validReads(0){};
	genbankTranscriptInfo(const std::string & initialiserText,std::string & transcript);
	std::string gene;
	int length,reads,validReads;
};

class genTranscriptInfoFile : public std::map<std::string,genbankTranscriptInfo> 
{
public:
	bool load(const std::string & filename);
};

#endif
