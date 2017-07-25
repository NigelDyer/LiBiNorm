// ***************************************************************************
// featureFile.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some functions for processing gtf and gff files
// ***************************************************************************

#ifndef FEATURE_FILE_H
#define FEATURE_FILE_H


#include <vector>
#include <string>
#include <map>
#include "containerEx.h"
#include "genbankFile.h"

extern bool htSeqCompatible;

class entryTag
{
public:
	char type1,type2;
	std::string val;
	entryTag(char type1,std::string val):type1(type1),type2(-1),val(val){};
	entryTag(char type1,char type2, std::string val) :type1(type1), type2(type2), val(val) {};
};


class entryTags : public std::vector<entryTag>
{
public:
	entryTags()/*:allowed(true)*/{};
	entryTags & operator = (const entryTags &tags)
	{
		std::vector<entryTag>::operator = (tags);
//		allowed = tags.allowed;
		return *this;
	};
	entryTags & operator += (const entryTags &tags)
	{
		for (const entryTag & i : tags)
			emplace_back(i);
//		if (!tags.allowed)
//			allowed = false;
		return *this;
	};
	entryTags(entryTags&&tags): std::vector<entryTag>(move(tags))/*, allowed(true)*/ { };

	static std::map<std::string,unsigned char> idMap ;
	static std::vector<std::string> revIdMap ;

	void add(const std::string type1,const std::string val);
	void add(const std::string type1, const std::string type2, const std::string val);

//	bool allowed;
};

class fileEntry;

typedef void(*parseTagsFunc)(fileEntry & fe, const char * start, size_t & len, const stringEx & type);
typedef bool(*printTagsFunc)(const entryTags & tags, outputDataFile * f);

//	This holds the information about how to parse the tages associated with the file format and the identity of any specific tag that is
//	being selected.  This is passed into the parser to be subsequently parsed by 
//		void parseval(const char * start, specificEntry & value,size_t & len) which is defined below
class specificEntry
{
public:
	specificEntry(fileEntry & fe, parseTagsFunc parseTags,const std::string type = "") :fe(fe),parseTags(parseTags), type(type) {};

	fileEntry & fe;
	parseTagsFunc parseTags;
	std::string type;
};

class qualType
{
public:
	qualType(int val=-1) :val(val) {};
	int val;
};


class fileEntry
{
public:

	fileEntry(fileEntry && entry):source(move(entry.source)), name(move(entry.name)), 
			id(move(entry.id)), chromosome(move(entry.chromosome)),
			type(move(entry.type)),	biotype(move(entry.biotype)), 
			strand(entry.strand),dot(entry.dot),
			valid(entry.valid),
			start(entry.start),finish(entry.finish),qual(entry.qual),
			tags(move(entry.tags)) {};
	fileEntry(const std::string & line, parseTagsFunc parseTags, const std::string & id_attribute = "");
	fileEntry(const std::string & chromosome,const std::string & featType,const featureSection & feature);
	~fileEntry() {};

	void addOffset(size_t begin, size_t end, int newLen);
	void output(printTagsFunc printTags,TsvFile & file);

	stringEx source,name;
	std::string id, chromosome, type;
	stringEx biotype;
	char strand,dot;
	bool valid;
	size_t start, finish;
	qualType qual;
	entryTags tags;

};



//	When we come to print a file entry then it is wrapped in this class so that the 
class fileEntryForPrinting
{
public:
	fileEntryForPrinting(const fileEntry & fe, printTagsFunc printTags) : fe(fe), printTags(printTags) {};

	const fileEntry & fe;
	printTagsFunc printTags;
};

class featureFile
{
	setEx<std::string> ignoredTranscriptTypes;
protected:
	//	For storing the headers of the gtf file as simple strings
	std::vector<std::string> headers;

	//	an entry for each line in the featureFile, indexed by chromosome name and chromosome position
	class featureMap : public  std::multimap<size_t, fileEntry>
	{
	public:
		ADD_ITER(position, feature);
	};

	class entryMapClass : public std::map<std::string, featureMap>
	{
	public:
		ADD_ITER(name, data);
		ADD_MAPPAIR(name, data);
	}entryMap;

	//	gff files index the entries by the chromsome ncbi identifier, so we have to
	//	use the chromosome lengths to indetify the user friendly name
	std::map<size_t, std::string> chromosomeMap;

	std::string filename;

	//	The tag format is different in a gff and gtf file so these are set to the functions required to parse and print them.
	parseTagsFunc parseTags;
	printTagsFunc printTags;

public:

	void ignoreTranscriptTypes(setEx<std::string> types) { ignoredTranscriptTypes.add(types); };
	void addToChromosomeMap(size_t length, const std::string & name) { chromosomeMap.emplace(length, name); };

	const std::map<std::string, featureMap> & entries(){ return entryMap;};

	//	For adding a shift to the values in a gtf file, when modifying the gene annotation based on reads from a specific sample
	void addOffset(const std::string & chromosome,int begin,int end, int newLen) ;
	void addGene (const std::string & chromosome,geneInfo & I);

	bool open(const stringEx & filename,const std::string & id_attribute = "",
		const std::string & feature_type = "",bool mapFile = false);
	bool save(const std::string & filename);

	//	A R value method for adding a parsed gtfEntry to the map
	void addEntry (fileEntry && entry);

	bool printEntries(const std::string & filename);

};

inline bool printVal(outputDataFile * f, const qualType & qual)
{
	if (qual.val == -1)
		fputs(".",f->fout);
	else
		fprintf(f->fout, "%i", qual.val);
	return true;
};

inline bool printVal(outputDataFile * f, const fileEntryForPrinting & fep)
{
	f ->printStart(fep.fe.chromosome,fep.fe.source,fep.fe.type,fep.fe.start,
		fep.fe.finish,fep.fe.dot,fep.fe.strand, fep.fe.qual);
	f->printSep();
	fep.printTags(fep.fe.tags,f);
	return true;
};



//	Support functions used by the parser class for parsing the gtfEntry text
namespace parserInternal
{
	void parseval(const char * start, specificEntry & value,size_t & len);
	void parseval(const char * start, qualType & value, size_t & len);
}

#endif
