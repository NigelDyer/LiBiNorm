// ***************************************************************************
// printEx.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For printing text files
// ***************************************************************************

#ifndef PRINTEX_H
#define PRINTEX_H

#include <stdio.h>
#include <stdarg.h>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include "inQuotes.h"

/*	Classes for printing data as variables separated by a specific character, e.g. comma or tabs.  The format is

printTsv(f,val1,val2,val3...);

Each datatype to be printed requires a printVal method.  More can be added eg to print custom classes.
printTsv puts a cr at the end, printTsvPartial prints the passed values without a terminating line feed.
*/

//	The core class that just handles the basic file operations
class outputFile
{
public:
	outputFile() :fout(0) {};
	~outputFile();

	//	Opening using a blank filename causes the output to be directed to stdout
	bool open(const char * filename = nullptr);
	bool open(const std::string & filename);
	bool is_open() { return (fout != 0); };

	//	Called automatically be the destructor
	void close();

	int printf(const char * line, ...);
	void flush();

//	operator FILE *() { return fout; };

	FILE * fout;

};


//	The workhorse.
//	A line of text can be printed using
//	outputDataFile file(":");		A colon separated file
//	file.open("filename");
//	file.print(a,b,c);				Prints a line in one go
//	file.printStart(a);
//	file.printMiddle(b);
//	file.printEnd(c);
class outputDataFile : public outputFile
{
protected:
	outputDataFile(const char * sep) : separator(sep),needSeparator(false){};

	//	Used internally to print each item in the list, managing the printing of separators
	void printInternal() {};

	template<typename First,typename... Rest>
	void printInternal(const First & value,const Rest & ... rest)
	{
		printSep();
		needSeparator = printVal(this,value);
		printInternal(rest...);
	};
public:
	void printSep(bool print = true)
	{
		if (separator && needSeparator && print)
			fputs(separator,fout);
		needSeparator = true;
	}

	template<typename intType>
	void printRepeat(intType N,const char * s = "")
	{
		for (intType i = 0; i < N; i++)
			printInternal(s);
	}

	//	Use this if you want to print a part line starting from the beginning.
	template<typename... Rest>
	void printStart(const Rest & ... rest)
	{
		needSeparator = false;
		printInternal(rest...);
	};


	//	Use this if you want to print a part line starting and finishing in the middle of a line.  You then need an 
	//	explicit xx.print() for the line feed at the end.
	template<typename... Rest>
	void printMiddle(const Rest & ... rest)
	{
		printInternal(rest...);
	};
	template<typename... Rest>

	void printEnd(const Rest & ... rest)
	{
		printInternal(rest...);
		print();
	};

	void print()
	{
		needSeparator = false;
		fprintf(fout,"\n");
	}

template<typename... Rest>
	void print(const Rest & ... rest)
	{
		printStart(rest...);
		print();
	};

private:
	const char * separator;
	bool needSeparator;

};

//	Some standard custom versions such as a tab separated file
class TsvFile : public outputDataFile
{
public:
	TsvFile() : outputDataFile("\t"){};
};

//	A comma separated file
class CsvFile : public outputDataFile
{
public:
	CsvFile() : outputDataFile(","){};
};

//	And a file with no separator, ie a simple text file
class TextFile : public outputDataFile
{
public:
	TextFile() : outputDataFile(0){};
};


/*************************************************************************
And now printVals for printing each of the data types.  Additional ones can be added 
for additional data types, either here if it is generic, or after the class is declared
for specific classes.  
*/
inline bool printVal(outputDataFile * f,int value)
{
	fprintf(f->fout,"%i",value);
	return true;
};
inline bool printVal(outputDataFile * f, short value)
{
	fprintf(f->fout, "%i", value);
	return true;
};
inline bool printVal(outputDataFile * f, long value)
{
	fprintf(f->fout, "%li", value);
	return true;
};

inline bool printVal(outputDataFile * f,unsigned int value)
{
	fprintf(f->fout,"%u",value);
	return true;
};
inline bool printVal(outputDataFile * f, unsigned long value)
{
	fprintf(f->fout, "%lu", value);
	return true;
};
#ifdef _WIN64
//	This is the same as unsigned long with the gcc compiler
inline bool printVal(outputDataFile * f,size_t value)
{
	fprintf(f->fout,"%zu",value);
	return true;
};
#endif

inline bool printVal(outputDataFile * f,double value)
{
    if ((value != 0) && ((fabs(value)<0.01) || (fabs(value)>1000000.0)))
		fprintf(f->fout, "%e", value);
	else
		fprintf(f->fout,"%f",value);

	return true;
};

inline bool printVal(outputDataFile * f,const char * value)
{
	fputs(value,f->fout);
	return true;
};
inline bool printVal(outputDataFile * f,const char & value)
{
	fprintf(f->fout,"%c",value);
	if (value == '\n') return false;    //No separator after newline
	return true;
};
inline bool printVal(outputDataFile * f,const std::string & value)
{
	fputs(value.c_str(),f->fout);
	return true;
};
inline bool printVal(outputDataFile * f,bool value)
{
	fprintf(f->fout, value?"Y":"N");
	return true;
};

//	For printing out C style arrays of things
template <typename _T,size_t N>
inline bool printVal(outputDataFile * f,_T (&value)[N])
{
	bool sep = false;
	for (size_t i = 0; i < N; i++)
	{
		f ->printSep(sep);
		printVal(f,value[i]);
		sep = true;
	}
	return sep;
}



template <class _Q>
inline bool printVal(outputDataFile * f,const fmtPrivate<_Q> & value)
{
	fprintf(f->fout,value.format,value.value);
	return true;
};

/*	
Some composite printVals for vectors and sets
*/

template <template <typename,typename> class Container, typename _T>
inline bool printVal(outputDataFile * f,const Container<_T,std::allocator<_T> > & value)
{
	bool sep = false;
	for (auto & v : value)
	{
		f ->printSep(sep);
		printVal(f,v);
		sep = true;
	}
	return sep;
};

template <template <typename,typename> class Container, typename _T>
inline bool printVal(outputDataFile * f,const fmtPrivate<Container<_T,std::allocator<_T> > > & value)
{
	bool sep = false;
	for (auto & v : value.value)
	{
		f ->printSep(sep);
		printVal(f,fmt(value.format,v));
		sep = true;
	}
	return sep;
};



//	Use inQuotes(value) to put quotes around something.

template <class _T>
inline bool printVal(outputDataFile * f,const inQuotesPrivate<_T> & value)
{
	fprintf(f->fout,"\"");
	printVal(f,value);
	fprintf(f->fout,"\"");
	return true;
};


/*****************************************************************************************
The default is that integers that are zero are printed.   If they are not to be printed then blankZero or _Z should be used, 
e.g. file.print(blankZero(X),Y); Again, using a template factory function to create the class means that the type is implicit and taken from the methods
parameters in a way that cannot be done with the class constructor.

TsvFile f;
f.open("filename");
f.print(val1,blankZero(val2),val3...);

*/

template <class _Q>
class blankZeroPrivate
{
//	Oddly, template<> needs to come before friend
template <class _R> friend blankZeroPrivate<_R> blankZero(const _R & val);
template <class _R> friend blankZeroPrivate<_R> _Z(const _R & val);
blankZeroPrivate(_Q v):value(v){};
public:
	_Q value;
};

template <class _R>
inline	blankZeroPrivate<_R> blankZero(const _R & val)
{
	return blankZeroPrivate<_R>(val);
};
template <class _R>
inline	blankZeroPrivate<_R> _Z(const _R & val)
{
	return blankZeroPrivate<_R>(val);
};


inline bool printVal(outputDataFile * f, const blankZeroPrivate<short> & value)
{
	if (value.value != 0)
		fprintf(f->fout, "%i", value.value);
	return true;
};
inline bool printVal(outputDataFile * f, const blankZeroPrivate<long> & value)
{
	if (value.value != 0)
		fprintf(f->fout, "%li", value.value);
	return true;
};

inline bool printVal(outputDataFile * f, const blankZeroPrivate<int> & value)
{
	if (value.value != 0)
		fprintf(f->fout, "%i", value.value);
	return true;
};

#ifdef _WIN64
inline bool printVal(outputDataFile * f,const blankZeroPrivate<size_t> & value)
{
	if (value.value != 0)
		fprintf(f->fout,"%zu",value.value);
	return true;
};
#endif
inline bool printVal(outputDataFile * f,const blankZeroPrivate<unsigned int> & value)
{
	if (value.value != 0)
		fprintf(f->fout,"%u",value.value);
	return true;
};

inline bool printVal(outputDataFile * f,const blankZeroPrivate<double> & value)
{
	if (value.value != 0)
		fprintf(f->fout,"%f",value.value);
	return true;
};
//	For bools only print Y for true, blank for false
inline bool printVal(outputDataFile * f, const blankZeroPrivate<bool> & value)
{
	if (value.value)
		fprintf(f->fout, "Y");
	return true;
};

//	This allows a blankZero to prefix a vector of things
template <template <typename,typename> class Container, typename _T>
inline bool printVal(outputDataFile * f,const blankZeroPrivate<Container<_T,std::allocator<_T> > > & value)
{
	bool sep = false;
	for (auto & v : value.value)
	{
		f ->printSep(sep);
		printVal(f,blankZero(v));
		sep = true;
	}
	return sep;
};



#endif
