// ***************************************************************************
// genbankFile.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Code for processing genbank files
// ***************************************************************************


#include "libCommon.h"
//	GenbankFile.h needs to be before parser.h so that the qualInt parseval can be seen and used by the gcc compiler.
//	The Microsoft compiler can work out the connection in either order
#include "genbankFile.h"
#include "parser.h"
#include "stringEx.h"

using namespace std;

bool printVal(outputDataFile * f,const orfInfo & value)
{
	f ->printStart(value.score,value.start,value.finish,value.finish2);
	return true;
};


int getVariant(const string & description)
{
	size_t pos;
	if ((pos = description.find("transcript variant")) != string::npos)
	{
		pos += 18;
		pos = description.find_first_not_of(" ",pos);
		if (description[pos] == 'X')
			pos++;
		return atoi(description.c_str() + pos);
	}
	else
		return 99;
}


string qualInt::asString()
{
	char res[100];
	if (qualifier != ' ')
		sprintf(res,"%c%i",qualifier,value);
	else
		sprintf(res,"%i",value);
	return res;
}
/*
extendableQuotedString & extendableQuotedString::operator = (const std::string & _S)
{
	_V.operator = (_S);
	return *this;
};
extendableQuotedString & extendableQuotedString::operator += (const std::string & _S)
{
	_V.operator += (_S);
	return *this;
};
*/

void parserInternal::parseval(const char *& start,qualInt & value,size_t & len)
{
	switch (*start)
	{
	case '>': case '<': 
		value.qualifier = *start;
		start++;len--;
		break;
	default:
		value.qualifier = ' ';
	}
	parserInternal::parseval(start,value.value,len);
}
void parserInternal::parseval(const char * start,extendableRemainderOfLine & value,size_t & len)
{
	value += start;
	len = value.size();
}

void parserInternal::parseval(const char *& start,extendableQuotedString & value,size_t & len)
{
	//	retlen is the length that is returned to ensure that the calling code skips
	// to just past the end of the string that has been found
	size_t skipQuotes = 0;

	if ((*start) == '\"')
	{
		//	If its a quoted string then look for an end quote and modify the 
		//	length accordingly if it is found
		start++;
		len--;
		const char * strEnd = strchr(start,'\"');
		if (strEnd)
		{
			len = strEnd - start;
			skipQuotes = 2;
		}
	}
	else
	{
		//	No quote at the start.   Is there a quote before the end, in which case
		//	stop the string just before the quote
		const char * strEnd = strchr(start+1,'\"');
		if ((strEnd) && ((size_t)(strEnd - start) < len))
		{
			len = strEnd - start;
			skipQuotes = 1;
		}
	}
	value += string(start,len);
	len += skipQuotes;
}

#define MAX_LOOKBACK 20

bool nucleotideSequence::seqMatch(size_t offset,const char * seq)
{
	for (size_t i = 0;i < strlen(seq);i++)
		if (toupper(at(i+offset)) != toupper(seq[i]))
			return false;
	return true;
}

void nucleotideSequence::getSecondORFs(size_t orf1End,int geneEnd,vector<orfInfo> & orfs)
{
	//	For each of the other two possible frames
	for (int v = 0;v < 2;v++)
	{
		int pos = (int)orf1End -4 - v;
		//	Look for all of the start codons going back MAX_LOOKBACK = 20 codons
		for (int i = 0; (i < MAX_LOOKBACK) && (pos > 0) ;i++,pos -=3)
		{
			if (seqMatch(pos,"atg"))
			{
				orfs[v].score += ((double)(MAX_LOOKBACK-i)/MAX_LOOKBACK);
				//	register the first start position
				if (orfs[v].start == 0)
					orfs[v].start = pos+1;	//positions are 1 based.  RNA sequence index is 0 based
			}
		}
		
		if(orfs[v].start)
		{
			//Look for the stop codon
			for(pos = orfs[v].start +2;pos < (int)size()-3;pos += 3)
			{
				if (seqMatch(pos,"t"))
				{
					if ((seqMatch(pos+1,"aa")) || (seqMatch(pos+1,"ag")) || (seqMatch(pos+1,"ga")))
					{
						orfs[v].finish = pos+3;
						break;
					}
				}
				if (pos > geneEnd)
				{
					orfs[v].finish = pos+3;
					break;
				}
			}
			//	And then look for a second stop codon
			pos += 3;
			for(;(pos < (int)size()-3 ) && (pos > geneEnd);pos += 3)
			{
				if (This[pos] == 't')
				{
					if ((seqMatch(pos+1,"aa")) || (seqMatch(pos+1,"ag")) || (seqMatch(pos+1,"ga")))
					{
						orfs[v].finish2 = pos+3;
						break;
					}
				}
			}

		}
		if (orfs[v].finish  < (int)orf1End)
		{
			orfs[v] = orfInfo();
		};

	}
}

size_t nucleotideSequence::getSecondStopCodon(size_t orf1Start)
{
	for(size_t pos = orf1Start+3;pos < min(size()-4,orf1Start+100);pos += 3)
	{
		if (seqMatch(pos,"t"))
		{
			if ((seqMatch(pos+1,"aa")) || (seqMatch(pos+1,"ag")) || (seqMatch(pos+1,"ga")))
			{
				return pos;
			}
		}
	}
	return 0;
}

featureSection::featureSection(genbankFile * f,const string & type,const char *& c) :type(type),linkedCDS(-1),variant(0)
{
	//to do...   This code does not cope with the join(1..2,3..5) syntax yet!
	parser(c," .",start,finish);

	c = f ->getLine();
	while ((c) && (c[5] == ' '))
	{
		string type;
		extendableQuotedString value;
		parser(c," /=\n",type,value);

		while ((c = f -> getLine()))
		{
			if ((c[5] != ' ') || (c[21] == '/'))
				break;
			parser(c," /=\n",value);
		}
		if (type == "gene")
		{
			size_t b = value.find('(');
			if (b != string::npos)
			{
				while (value[b-1] == ' ')
					b--;
				value = value.substr(0,b);
			}
		}

		emplace(type,value);
	}
}

const string & geneInfo::getVal(const char * str1,const char * str2)
{
	static string nullString;
	for (featureSection & f : features)
	{
		if (f.type == str1)
		{
			auto q = f.find(str2);
			if (q == f.end())
				return nullString;
			return q->second;
		}
	}
	return nullString;
}

const string & geneInfo::ncRNA_class()
{
	return getVal("ncRNA","ncRNA_class");
}
const char * geneInfo::geneName()
{
	return getVal("gene","gene").c_str();
}


const char * geneInfo::readNext(genbankFile * f)
{
	features.clear();
	type.clear();
	string featType;

	const char * c = f -> line();

	if (c[0] != ' ')
		return 0;
	do
	{
		parser(c," \n",featType);
		if (featType == "gene")
		{
			features.push_back(featureSection(f,featType,c));
			break;
		}
		else
			featureSection fs(f,featType,c);
	} while (c);

	if (c == nullptr)
		return c;

	if (strncmp(c,"ORIGIN",6) == 0)
		return c;

	parser(c," \n",featType);

	do
	{
		static setEx<string> supportedFeatures {{"CDS","ncRNA","mRNA"}};
		if (supportedFeatures.contains(featType))
		{
			features.push_back(featureSection(f,featType,c));
			if (featType == "CDS")
			{
				for (featureSection & f: features)
				{
					if ((f.type == "mRNA") && (f.linkedCDS == -1))
						f.linkedCDS = (int)features.size()-1;
				}
			}
		}
		else
			featureSection fs(f,featType,c);
		if (type.empty()) 
			type = featType;


		if ((c) && (c[0] != ' '))
			return c;
		parser(c," \n",featType);
	} while ((c != 0) && (featType != "gene"));

	return c;
}

const string & genbankLocus::mol_type()
{
	return source["mol_type"];
}

bool genbankLocus::readNext(genbankFile * f)
{
	LOCUS_length = variant = 0;
	VERSION = LOCUS_id = DEFINITION = SOURCE = "";
	source.clear();
	nucSeq.clear();
//static bool found = false;
	const char * c = f -> getLine();
	if (c == nullptr)
		return false;


	do {
		string entryType;
		parser(c," \n",entryType);

		if (entryType == "LOCUS") {
			parser(c," \n",LOCUS_id,LOCUS_length);
			c = f -> getLine();
		}
		else if (entryType == "VERSION") {
			vector<string> params;
			parser(c," :\n",VERSION,params);
			for (size_t i = 0;i < params.size()-1;i++)
				if (params[i] == "GI") GI = params[i+1];
			c = f -> getLine();
		}
		else if (entryType == "DEFINITION") {
			parser(c," \n",DEFINITION);
			c = f -> getLine();
			while ((c) && (c[0] == ' '))
			{
				//	Deal with multiline definitions
				parser(c," \n",DEFINITION);
				c = f -> getLine();
			}
			//	Extract the transcript variant info, which is after the text transcript variant and will be either <n> or X<n>

			variant = getVariant(DEFINITION);
			DEFINITION.truncateFrom("(");

		}
		else if (entryType == "SOURCE") {
			parser(c," \n",SOURCE);
			SOURCE.truncateFrom("(");
				
			if ( strncasecmp(SOURCE.c_str(),DEFINITION.c_str(),SOURCE.size())== 0)
				DEFINITION = DEFINITION.substr(SOURCE.size());

			DEFINITION[0] = toupper(DEFINITION[0]);
		}
		else if (entryType == "FEATURES") {
			c = f ->getLine();
			string featType;
			parser(c," \n",featType);
			if (featType == "source")
				source = featureSection(f,featType,c);
			return true;
		}
		else
			c = f -> getLine();

	} while (c != nullptr);


	return (c != 0);
}

bool genbankLocus::readOriginSeq(genbankFile * f)
{
	const char * c = f -> line();

	do {
		string entryType;
		parser(c," \n",entryType);
		if (entryType == "ORIGIN") {
			while((c = f -> getLine()) != 0)
			{
				if (strncmp(c,"//",2) == 0)
					return true;
				int pos;
				vector<string> seqs;
				parser(c," \n",pos,seqs);
				for (string & seq : seqs)
					nucSeq += seq;
			}
		}
	} while (true);
	return true;
}



bool  genbankFile::getNextLocus(genbankLocus & locus)
{
	return locus.readNext(this);
}

bool  genbankFile::getNextGene(geneInfo & gene)
{
	return (gene.readNext(this) != 0);
}


bool genbankFile::open(const char * filename)
{
	f.open(filename);
	return (f.is_open());
};


/**************************************************************
	For processing sets of genBank data as retreived from a genbank file

*/

genbankInfoParser::genbankInfoParser(const stringEx & header)
{
	vector<string> colHeaders;
	parseTsv(header,colHeaders);
	columns.resize(colHeaders.size(),null);
	for (size_t i = 0;i < colHeaders.size();i++)
	{
		if (colHeaders[i] == "Ensembl Gene ID")
			columns[i] = geneIdPos;
		else if (colHeaders[i] == "Associated Gene Name")
			columns[i] = geneNamePos;
		else if (colHeaders[i] == "Gene Start (bp)")
			columns[i] = startPos;
		else if (colHeaders[i] == "Gene End (bp)")
			columns[i] = endPos;
		else if (colHeaders[i] == "Chromosome Name")
			columns[i] = chromosomePos;
		else if ((colHeaders[i] == "Gene type") || (colHeaders[i] == "Gene Biotype"))
			columns[i] = typePos;
	}
}

genbankInfo::genbankInfo(const char * initialiserString)
{
	parseTsv(initialiserString,type,gene,definition,variant,ncRNAtype,length,chromosome,geneStart,CDSstart,CDSfinish,finish2,
			orf2offset,orf2score,orf2start,orf2finish,orf2finish2,orf2length,reads,validReads,variants);
}


genbankInfo::genbankInfo(const std::string & initialiserString,const genbankInfoParser & info,std::string & id)
{
	vector<string> values;
	parseTsv(initialiserString,values);
	for (size_t i = 0; i < min(info.columns.size(),values.size());i++)
	switch (info.columns[i])
	{
	case genbankInfoParser::geneIdPos: id = values[i]; break;
	case genbankInfoParser::geneNamePos: gene = values[i];break;
	case genbankInfoParser::typePos: type = values[i];break;
	case genbankInfoParser::chromosomePos: chromosome = values[i];break;
	case genbankInfoParser::startPos: geneStart = atoi(values[i].c_str());break;
	case genbankInfoParser::endPos: length = atoi(values[i].c_str()) - geneStart;break;
	default:;
	}
}
/*
genbankTranscriptInfo::genbankTranscriptInfo(const string & initialiserText,string & transcript)
{
	parseTsv(initialiserText,transcript,gene,length,reads,validReads);
}
*/

bool genInfoFile::load(const string & filename)
{
	ifstream m_file;
	m_file.open(filename);
	if (!m_file.is_open()) 
		return false;

	string line;
	std::getline(m_file,line);
	while (	std::getline(m_file,line))
	{
		string root,version;
		const char * c = line.c_str();
		//When passing in a const char * to parseTsv it moves the pointer to the end of the last parameter that was read
		parseTsv(c,root,version);
		emplace(version,genbankInfo(c));
	}
	for (auto gbi : This)
		chromosomeInfo[gbi.second.chromosome].emplace(gbi.second.CDSstart.value,gbi.first);

	return true;

}
/*
bool genTranscriptInfoFile::load(const string & filename)
{
	ifstream m_file;
	m_file.open(filename);
	if (!m_file.is_open()) 
		return false;

	string line;
	std::getline(m_file,line);	//For title
	while (	std::getline(m_file,line))
	{
		string transcript;
		genbankTranscriptInfo gbti(line,transcript);
		emplace(transcript,gbti);
	}

	return true;

}
*/
