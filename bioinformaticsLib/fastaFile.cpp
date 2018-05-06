// ***************************************************************************
// fastaFile.cpp (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 28 February 2018
// ---------------------------------------------------------------------------
// Some functions for processing fasta files
// ***************************************************************************

#include <math.h>
#include <algorithm> 
#include "printEx.h"
#include "libParser.h"
#include "fastaFile.h"

using namespace std;
/*


  #           (bases used for these tables, for reference)
   # 1 TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
   # 2 TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
   # 3 TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
*/

sequence sequence::toAA() const
{
	static const char * aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	sequence AA;
	int index = 0;
	int nCount = 0;
	bool X = false;
	for (auto c : This)
	{
		switch (c)
		{
		case 'T': case 't': index *= 4;break;
		case 'C': case 'c': (index *= 4) += 1;break;
		case 'A': case 'a': (index *= 4) += 2;break;
		case 'G': case 'g': (index *= 4) += 3;break;
		default: X = true;
		}
		if (++nCount == 3)
		{
			if (X)
			{
				AA += 'X';
				X = false;
			}
			else
				AA += aa[index];
			nCount = 0;
			index = 0;
		}
	}
	return AA;
};


sequence sequence::complement() const
{
	size_t len = size();
	sequence s(len,'X');
	s.comp = !comp;
	for (size_t i = 0;i < size();i++)
	{
		switch (at(i))
		{
		case 'A': s[len-i-1] = 'T';break;
		case 'C': s[len-i-1] = 'G';break;
		case 'G': s[len-i-1] = 'C';break;
		case 'T': s[len-i-1] = 'A';break;
		case 'a': s[len-i-1] = 't';break;
		case 'c': s[len-i-1] = 'g';break;
		case 'g': s[len-i-1] = 'c';break;
		case 't': s[len-i-1] = 'a';break;
		case 'N': s[len-i-1] = 'N';break;
		case 'n': s[len-i-1] = 'n';break;
		}
	}
	return s;
}

sequence sequence::reversed() const
{
	size_t len = size();
	sequence s(len,' ');
	for (size_t i = 0;i < size();i++)
	{
		s[len-i-1] = at(i);
	}
	return s;
}


void trimR(string & str)
{
	if (str[str.length()-1] == 13)
		str.resize(str.length()-1);
}


template <> bool fastaFile<ofstream>::writeEntry(const string & sequence,const char * header)
{
	m_file << ">" << header << endl;
	const char * cursor  = sequence.c_str();
	const char * end = cursor + sequence.length();

	while (cursor < end)
	{
		ptrdiff_t len = min<ptrdiff_t>(LINE_LEN,end-cursor);
		m_file.write(cursor,len);
		m_file << endl;
		cursor += LINE_LEN;
	}
	sections.push_back(fastaSection(header,(int)sequence.size()));
	return true;
}

template <> void fastaFile<ifstream>::close(bool index)
{
	if (m_file.is_open())
		m_file.close();
}
template <> void fastaFile<ifstream>::getline(string & line)
{
	std::getline(m_file,line);
}
template <> void fastaFile<ofstream>::close(bool index)
{
	if (m_file.is_open())
	{
		m_file.close();
		if (index)
		{
			int cumulative_length = 0;
			TsvFile f;
			f.open(filename + ".fai");
			m_file.open((filename + ".fai").c_str(),	ios_base::binary);
			for (auto & section : sections)
			{
				int lineLen = min(LINE_LEN,section.length);
				f.print(section.name,section.length,cumulative_length+section.name.size()+3,lineLen,lineLen + 2);
				cumulative_length += (section.length + ((int)ceil(1.0 * section.length-1)/LINE_LEN) * 2 + (int)section.name.size()+5);
			}
			f.close();

/*			m_file.open((filename + ".fai").c_str(),	ios_base::binary);
			for (auto & section : sections)
			{
				m_file << section.name << "\t" << section.length << "\t" << cumulative_length + section.name.size()+3 << 
					"\t" << LINE_LEN << "\t" << (LINE_LEN + 2) << "\n";
				cumulative_length += (section.length + (int)ceil((1.0 * section.length)/LINE_LEN) + (int)section.name.size()+1);
			}
			m_file.close();*/
		}
	}
}


bool fastaFileRead::readName()
{
	if (m_file.eof())
		return false;
	char c = m_file.get();
	if (c != '>') 
		return false;

	getline(m_name);
//	trimR(m_name);
	return true;
}

bool fastaFileRead::readEntry()
{
	if (!readName())
		return false;

	getline(m_sequence);
//	trimR(m_sequence);

	return true;
}

bool fastaFileRead::readEntry(string & sequence,bool minLoad)
{
	if (!readName())
		return false;

	sequence.resize(0);
	while((m_file.peek() != '>') && (!m_file.eof()))
	{
		std::getline(m_file,m_sequence);
		if(!minLoad)
		{
			sequence += m_sequence;
		}
	}
	return true;
}

bool fastaFileRead::readAll()
{
	string sequence;
	while (readEntry(sequence))
	{
		sequences.emplace(m_name.substr(0,m_name.find_first_of(' ')),sequence);
	}
	return true;
}

const string & fastaFileRead::getSeq(const string & id)
{
	auto sequenceIterator = sequences.find(id);
	if (sequenceIterator != sequences.end())
		return sequenceIterator -> second;
	static const string dummy;
	return dummy;
}

bool fastqFile::getNext(string & seq,string & name)
{
	if (m_file.eof())
		return false;
	string line;
	while((m_file.peek() != '@') && (!m_file.eof()))
	{
		getline(m_file,line);
	}
	getline(m_file,line);
	if (line.size() == 0)
		return false;
	size_t pos = line.find(' ');
	if (pos != string::npos)
		name = line.substr(1,pos-1);
	else
		name = line.substr(1);
	getline(m_file,seq);
	getline(m_file,line);
	getline(m_file,line);
	return true;

}




