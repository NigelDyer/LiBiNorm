// ***************************************************************************
// fastaFile.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some functions for processing fasta files
// ***************************************************************************

#ifndef FASTA_FILE_H
#define FASTA_FILE_H

#define MAX_LINE 5000
#define LINE_LEN 80

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <memory>


void trimR(char * buffer);

struct sequence : public std::string
{
protected:
	bool comp;
public:
	sequence(const char * s):std::string(s),comp(false) {};
	sequence():comp(false){};
	sequence(size_t len, char c):std::string(len,c),comp(false){};
	sequence(const std::string & s):std::string(s),comp(false){};
	sequence(const sequence & s):std::string(s),comp(s.comp){};
	sequence complement() const;
	sequence reversed() const;
	sequence toAA() const;
};




//	Can be either an ifstream or an ofstream, depending on whether the fastafile is being written or read.
template <class stream> class fastaFile 
{
public:
	enum 
	{
		FORWARD = true,
		REVERSE = false
	};
protected:
	stream m_file;
	int max_chromosome;
	std::string filename;

	std::string m_name;
	std::string m_sequence;

	class fastaSection
	{
	public:
		fastaSection(std::string s,int l):name(s),length(l){};
		std::string name;
		int length;
	};
	std::vector<fastaSection> sections;
public:
	fastaFile();
	~fastaFile();
	const char * Name() {return m_name.c_str();};
	std::string & NameStr() { return m_name;};
	bool eof() {return m_file.eof();};
	void getline(std::string & line);

	bool open(const char * filename,const char * mode = "");
	bool open(const std::string & filename,const char * mode = "");
	bool open(const std::ostringstream & filename,const char * mode = "");
	void close(bool index = false);
	bool writeEntry(const std::string & sequence,const char * header);
	bool writeEntry(const std::string & sequence,const std::string & header){ return writeEntry(sequence,header.c_str());} ;
};

template <class stream> fastaFile<stream>::fastaFile(void) :max_chromosome(0) {};

template <class stream> fastaFile<stream>::~fastaFile(void)
{
	close();
}

template <class stream> bool fastaFile<stream>::open(const char * name,const char * mode)
{
	filename = name;
//	ios_base::openmode m = ios_base::binary;
	m_file.open(name);
	return m_file.is_open();
}

template <class stream> bool fastaFile<stream>::open(const std::string & name,const char * mode)
{
	return open(name.c_str(),mode);
}

template <class stream> bool fastaFile<stream>::open(const std::ostringstream & name,const char * mode)
{
	return open(name.str(),mode);
}

class fastaFileRead: public fastaFile<std::ifstream>
{
public:
	std::map<std::string,std::string> sequences;
	fastaFileRead(){};
	bool readName();
	bool readEntry();
	bool readEntry(std::string & sequence,bool minLoad = false);
	bool readAll();
	const std::string & getSeq(const std::string & id); 
};

typedef fastaFile<std::ofstream> fastaFileWrite;

//	A fastafile that can be pushed onto vectors and inserted into maps  Use a typedef rather than a new class because
//	a new class would require rValue copy constructors etc.
typedef  std::shared_ptr<fastaFileRead> fastaFileReadEx;


class fastqFile
{
	std::ifstream m_file;
public:
	bool open(const std::string & filename)
	{
		m_file.open(filename.c_str(),std::ios_base::in);
		return m_file.is_open();
	}
	bool getNext(std::string & seq,std::string & name);

};

#endif