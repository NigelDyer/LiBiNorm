// ***************************************************************************
// parser.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 28 February 2018
// ---------------------------------------------------------------------------
// For parsing text files
// ***************************************************************************

#ifndef PARSER_H
#define PARSER_H

#include "libCommon.h"
#include <string>
#include <map>
#include <vector>
#include <fstream>

/*
	parseval method parses a sourceString into separate tokens
	
	parseval(sourceString,delimiter,value1,value2(,delim2));

	sourceString can either be a const char * or a std::string

	Currently supported value types are int, double,char, string.  The last paremeter can be vector of things in which case
	entries are added until the string is finished
*/
/*
	parsevals are put into a namespace to keep them separate.  I would have liked to have included them in the class
	The objective is to allow additional definitions to be added by users of the parser class.  I was hoping to do this
	by using inherited classes, but unfortunately template functions in child classes are not able to see methods in the 
	parent class, so any that the child wants to use have to be explicitly redeclared and the code then calls the parent method.
	However declaring additional methods in the child include file before including the parent definition does seem to allow the parent
	non-member methods to be compiled using methods defined in the child header
*/

class allowNulls
{
	friend class parser;
	const char * delimiter;
public:
	allowNulls(const char * delim) : delimiter(delim){};

};
namespace parserInternal
{
	//	Allows void* to be passed in to skip an unwanted value
	inline void parseval(const char * start,void * & value,size_t & len){};

	void parseval(const char * start,short & value, size_t & len);
	void parseval(const char * start,long & value, size_t & len);
	void parseval(const char * start,size_t & value, size_t & len);
	void parseval(const char * start,int & value,size_t & len);
	void parseval(const char * start,double & value,size_t & len);
	void parseval(const char * start,float & value,size_t & len);
	void parseval(const char * start,long double & value,size_t & len);
	void parseval(const char * start,char & value,size_t & len);
	void parseval(const char * start,std::string & value,size_t & len);
}
using namespace parserInternal;

class parser
{
	const char * delim;
	bool multiSkip;

protected:
	parser (const char * delimiter,bool multiSkip=false):delim(delimiter),multiSkip(multiSkip){}; 

	// A dummy method which is called when there are no further columns to process
	void parseNextEntry(const char * start) {};

	// The function itself, which is called recursively, each time peeling off and processing one entry 
	//  from the tab delimited line.  Once the string is exhausted the start pointer is left on the terminating null and this 
	//  calls the parseval method with the length set to zero, simulating an empty column
	//	This version steps over multiple instances of consecutive delimiters
	template <typename First,typename... Rest>
	void parseNextEntry (const char *& start,First & first,Rest& ... rest)
	{
		size_t len = 0;
		if (start)
		{
			if (multiSkip)
				start = start + strspn(start,delim);

			len = strcspn(start,delim);

			//  Look for the next delimiter.
			parseval(start,first,len);

			//  move the cursor, but don't move past the end of the string
			start += len;
			if (*start)
				start++;
		}
		else
			parseval(start,first,len);

		parseNextEntry(start,rest...);
	}

	template<class _Kty,class _Ty>
	void parseNextEntry (const char *& start,std::map<_Kty,_Ty> & data)
	{
		size_t len = 0;

		if (multiSkip)
			start = start + strspn(start,delim);
		len = strcspn(start,delim);

		//  Look for the next delimiter.
		int entries;
		parseval(start,entries,len);

		start += len;
		if (*start)
			start++;

		for (int i = 0;i < entries && (*start);i++)
		{
			if (multiSkip)
				start = start + strspn(start,delim);
			len = strcspn(start,delim);

			//  Look for the next delimiter.
			_Kty var;
			parseval(start,var,len);

			//  move the curser
			start += len;
			if (*start)
				start++;

			_Ty var2;
			parseval(start,var2,len);
			//	Provides misleading information in the case of a set, as 
			//	it suggest that it should search from the end to find the
			//	right location.  This does however allow common code t be used
			//	for sets and vectors.
			data.emplace(var,var2);
			//  move the cursor
//			start += len;
//			if (*start)
//				start++;
		}
	};


	//	Works for sets and vectors, puts everything else that is found into a set, vector setEx or vectorEx
	template <template <typename,typename> class Container, typename _T>
	void parseNextEntry (const char *& start,Container<_T,std::allocator<_T> > & first)
	{
		size_t len = 0;
		while (*start)
		{
			if (multiSkip)
				start = start + strspn(start,delim);
			len = strcspn(start,delim);

			//  Look for the next delimiter.
			_T var;
			parseval(start,var,len);

			//	Provides misleading information in the case of a set, as 
			//	it suggest that it should search from the end to find the
			//	right location.  This does however allow common code t be used
			//	for sets and vectors.
			first.insert(first.end(),var);
			//  move the cursor
			start += len;
			if (*start)
				start++;
		}

	}

	//	This is for C arrays of a set size e.g. std::string[4]
	template <typename _T,typename... Rest,size_t N>
	void parseNextEntry (const char *& start,_T (&first)[N],Rest&...rest)
	{
		for (size_t i = 0;i < N;i++)
		{
			size_t len = 0;
			if (*start)
			{
				if (multiSkip)
					start = start + strspn(start,delim);
				len = strcspn(start,delim);

				//  Look for the next delimiter.
				parseval(start,first[i],len);

				//	Provides misleading information in the case of a set, as 
				//	it suggest that it should search from the end to find the
				//	right location.  This does however allow common code t be used
				//	for sets and vectors.
				//  move the cursor
				start += len;
				if (*start)
					start++;
			}
		}
		parseNextEntry(start,rest...);

	}


	//
	//	A local delimiter can be passed in as a character string 
	template <typename Value>
	void parseNextEntry (const char *& start,Value & value,const char * localDelim)
	{
		size_t len = 0;
		if (start)
		{
			if (multiSkip)
				len = strcspn(start,localDelim);

			//  Look for the next delimiter.
			parseval(start,value,len);

			//  move the cursor
			start += len;
			start += strspn(start,localDelim);
		}
		else
			parseval(start,value,len);
	}

	template <typename First,typename... Rest>
	void parseNextEntry (const char *& start,First & first,const char * localDelim,Rest& ... rest)
	{
		size_t len = 0;
		if (start)
		{
			len = strcspn(start,localDelim);

			//  Look for the next delimiter.
			parseval(start,first,len);

			//  move the cursor
			start += len;
			start += strspn(start,localDelim);
		}
		else
			parseval(start,first,len);

		parseNextEntry(start,rest...);
	}
/*
	void parseNextEntry (const char *& start,vector<string> & rest)
	{
		string name;
		do
		{
			parseNextEntry(start,name);
			if(name.size())
				rest.push_back(name);
			else
				break;
		} while (true);
	}
*/
//	template <typename... Values>
//	parser (const char * delimiter,Values& ... values) : delim (delimiter),multiSkip(false){};

public:
	template <typename... Values>
	parser (const char *& start,const char * delimiter,Values& ... values) : delim (delimiter),multiSkip(true)
	{
		parseNextEntry(start,values...);
	}

	template <typename... Values>
	parser (const std::string & start,const char * delimiter,Values& ... values) : delim (delimiter),multiSkip(true)
	{
		const char * cursor = start.c_str();
		parseNextEntry(cursor,values...);
	}
	template <typename... Values>
	parser (const std::string & start,const allowNulls & delimiter,Values& ... values) : delim (delimiter.delimiter),multiSkip(false)
	{
		const char * cursor = start.c_str();
		parseNextEntry(cursor,values...);
	}

	
	template <typename... Values>
	parser (const std::string & start,Values& ... values) : delim (""),multiSkip(true)
	{
		const char * cursor = start.c_str();
		parseNextEntry(cursor,values...);
	}
};


class parseTsv : public parser
{
protected:
	parseTsv (): parser("\r\n\t"){};

public:
	template <typename... Values>
	parseTsv (const std::string & start,Values& ... values): parser("\r\n\t")
	{
		const char * cursor = start.c_str();
		parseNextEntry(cursor,values...);
	};
	template <typename... Values>
	parseTsv (const char *& start,Values& ... values): parser("\r\n\t")
	{
		parseNextEntry(start,values...);
	};
};

class parseCsv : public parser
{
protected:
	parseCsv (): parser("\r\n,"){};

public:
	template <typename... Values>
	parseCsv (const std::string & start,Values& ... values): parser("\r\n,")
	{
		const char * cursor = start.c_str();
		parseNextEntry(cursor,values...);
	};
	template <typename... Values>
	parseCsv (const char *& start,Values& ... values): parser("\r\n,")
	{
		parseNextEntry(start,values...);
	};
};


class parseTsvFile
{
public:
	bool open(const std::string & filename);
	void read() {
		if (!file.eof())
			std::getline(file, line);
	}
#ifdef LIBINORM
	template<typename Value, typename... Rest> void read(Value & value, Rest&... rest)
	{
		if (!file.eof())
		{
			std::getline(file, line);
			parseTsv(line, value);
			read(rest...);
		}
	}
#else
	template<typename... Values> void read(Values&... values)
	{
		if (!file.eof())
		{
			std::getline(file, line);
			parseTsv(line, values...);
		}
	}
#endif
	bool read(const char * expected)
	{
		if (file.eof())
			return false;
		{
			std::getline(file, line);
			if (line.substr(0, strlen(expected)) != expected)
				exitFail("Unexpected ", line);
		}
		return true;
	}

	template<typename Value1, typename Value2, typename... Rest> void read(std::map<Value1,Value2> & value, Rest&... rest)
	{
		while (!file.eof())
		{
			if (file.eof()) return;
			std::getline(file, line);
			if (line == "")
			{
				read(rest...);
				return;
			}
			Value1 value1;
			Value2 value2;
			parseTsv(line, value1,value2);
			value.emplace(value1, move(value2));
		}
	}
#ifdef LIBINORM
	template<typename Value, typename... Rest> void read(std::vector<Value> & value, Rest&... rest)
	{
		for (size_t i = 0;i < value.size();i++)
		{ 
			if (file.eof()) return;
			std::getline(file, line);
			parseTsv(line, value[i]);
		}
		read(rest...);
	}
#endif
	bool eof()
	{
		return file.eof();
	}

private:
	std::ifstream file;
	std::string line;
};




#endif
