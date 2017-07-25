// ***************************************************************************
// parser.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For parsing text files
// ***************************************************************************

#include "libCommon.h"
#include "parser.h"

using namespace std;

void parserInternal::parseval(const char * start,char & value,size_t & len)
{
	if (len)
		value = *start;
	else
		value = ' ';
}

void parserInternal::parseval(const char * start,int & value,size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try 
		{
			value = stoi(string(start,len));
		}
		catch(...)
		{
			value = 0;
		}
	}
}
void parserInternal::parseval(const char * start, long & value, size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try
		{
			value = stol(string(start, len));
		}
		catch (...)
		{
			value = 0;
		}
	}
}
void parserInternal::parseval(const char * start, size_t & value, size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try
		{
			value = stoul(string(start, len));
		}
		catch (...)
		{
			value = 0;
		}
	}
}
void parserInternal::parseval(const char * start, short & value, size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try
		{
			value = stoi(string(start, len));
		}
		catch (...)
		{
			value = 0;
		}
	}
}
void parserInternal::parseval(const char * start,float & value,size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try 
		{
			value = stof(string(start,len));
		}
		catch(...)
		{
			value = 0;
		}
	}
}
void parserInternal::parseval(const char * start,double & value,size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try 
		{
			value = stod(string(start,len));
		}
		catch(...)
		{
			value = 0;
		}
	}
}
void parserInternal::parseval(const char * start,long double & value,size_t & len)
{
	if (len == 0)
		value = 0;
	else
	{
		try 
		{
			value = stod(string(start,len));
		}
		catch(...)
		{
			value = 0;
		}
	}
}
/*
void parserInternal::parseval(const char * start,std::string & value,size_t & len)
{
	if (len == 0)
		value.clear();
	else
		value = string(start,len);
}
*/
void parserInternal::parseval(const char * start,string & value,size_t & len)
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
	value.assign(start,len);
	len += skipQuotes;
}



bool parseTsvFile::open(const std::string & filename)
{
	file.open(filename);
	return file.is_open();
}


