// ***************************************************************************
// stringEx.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// An extension to std::string that provides additional string handling functions
// ***************************************************************************

#ifndef STRINGEX_HEADER
#define STRINGEX_HEADER

#include <string.h>
#include <stdio.h>
#include <string>
#include <set>
#include <vector>
#include <sstream>
#include "inQuotes.h"
/*	A string extension class that supports a constructor that concatenates together all of the values passed to the constructor 
The stringEx class

*/
#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#define BUFFERSIZE  512
class stringEx : public std::string
{
public:
	stringEx (){};

	//	Define like this so that it is not used for the default constructor (which is above) because the Visual C++
	//	code checker does not spot a variadic constructor as a valid default contructor
	template <typename First,typename... Rest>
	stringEx (const First& first,const Rest& ... rest)
	{
		append(first,rest...);
	}

//  Additional types can be supported by adding an additional to_string operator.   See genomicPosition for an example.
//	Note that the to_string must be placed in the std:: namespace otherwise gcc will only see it and not look
//	for possibilities in the std namespace.   The microsoft compiler will give both equal weight so it does
//	not require to_string it to be in the gcc namespace
	template <typename _t>
		void append (const _t val) { std::string::append(std::to_string(val));};

private:

//	Dont define a null append function as this should unsure the compiler generates an error
//  if someone tries to append nothing, which is probbaly not what they intended to do
//	void append () {};
public:
	//	A few simple append types where we can use the std::string version
	void append (const stringEx & val) { std::string::append(val);};
	void append (const std::string & val) { std::string::append(val);};
	void append (const char * val) { std::string::append(val);};
	void append (char * const & val) { std::string::append(val);};
//	void append (const char & val) { push_back(val);};
	void append (const char val) { push_back(val);};

	operator bool() const { return size()>0;};

	template <typename _R>
	void append (const inQuotesPrivate<_R> & val) 
	{ 
		push_back('\"');
		append(val.value);
		push_back('\"');
	};

	template <typename _R>
	void append (const fmtPrivate<_R> & val) 
	{ 
		char buff[BUFFERSIZE];
		snprintf(buff,BUFFERSIZE,val.format,val.value);
		append(buff);
	};

	template <typename _t>
		stringEx & operator += (const _t & val)
		{
			 append(val);
			 return *this;
		}

	template <typename First,typename... Rest>
	stringEx & append (const First & first,const Rest& ... rest)
	{
		append(first);
		append(rest...);
		return *this;
	}

/*	template <typename First>
	stringEx & append (const First & first)
	{
		append(first);
		return *this;
	}
	*/
	stringEx & appendTsv (){ return *this;};

	template <typename First>
	stringEx & appendTsv (const First & first)
	{
		append(first);
		return *this;
	}


	template <typename... Rest>
	stringEx & appendTsv (const char first,const Rest& ... rest)
	{
		append(first);
		if (first != '\n')
			push_back('\t');
		appendTsv(rest...);
		return *this;
	}
	template <typename First,typename... Rest>
	stringEx & appendTsv (const First & first,const Rest& ... rest)
	{
		append(first);
		push_back('\t');
		appendTsv(rest...);
		return *this;
	}


	//	Works for vectors, sets and multisets of things
	template <template <typename,typename> class Container, typename _T>
	stringEx & appendTsv (const Container<_T,std::allocator<_T> > & Data)
	{
		bool firstParam = true;
		for (auto & _V : Data)
		{
			if (firstParam)
				firstParam = false;
			else
				push_back('\t');
			append(_V);
		}
		return *this;
	}

	//	Works for vectors, sets and multisets of things
	template <template <typename,typename> class Container, typename _T>
	stringEx & appendWithSep(const Container<_T,std::allocator<_T> > & Data,char sep)
	{
		bool firstParam = true;
		for (auto & _V : Data)
		{
			if (firstParam)
				firstParam = false;
			else
				push_back(sep);
			append(_V);
		}
		return *this;
	}

		//	Works for vectors, sets and multisets of things which are followed by other items
	template <template <typename,typename> class Container, typename _T,typename... Rest>
	stringEx & appendTsv (const Container<_T,std::allocator<_T> > & Data,const Rest& ... rest)
	{
		appendTsv(Data);
		push_back('\t');
		appendTsv(rest...);
		return *this;
	}

	template <typename... Rest>
	stringEx replaceSuffix(const Rest& ... rest) const
	{
		return stringEx(substr(0,find_last_of('.')),rest...);
	}
	template <typename... Rest>
	stringEx insertAtEndOfFileRoot(const Rest& ... rest) const
	{
		return stringEx(substr(0,find_last_of('.')),rest...) +  substr(find_last_of('.'));
	}

	stringEx removeSuffix() const
	{
		return substr(0,find_last_of('.'));
	}
	stringEx suffix() const
	{
		return substr(find_last_of('.')+1);
	}
	template <typename... Rest>
	stringEx replaceFilename(const Rest& ...  newName) const 
	{
		return stringEx(substr(0,find_last_of("/\\")+1),newName...);
	}
	stringEx filename() const
	{
		return substr(find_last_of("/\\")+1);
	}
	stringEx directory() const
	{
		return substr(0,find_last_of("/\\") + 1);
	}
	stringEx toLower() const;

	/*	stringEx & replaceAll(const std::string & from,char to)
	{
		size_t f = 0;
		while ((f = find_first_of(from,f)) != std::string::npos)
			at(f) = to;
		return *this;
	}	*/
	bool contains(const char * _V) const
	{
		return (find(_V) != npos);
	}
	bool endsWith(const std::string _V) const
	{
		return (substr(size()-_V.size()) == _V);
	}
	bool startsWith(const std::string _V) const
	{
		return (substr(0,_V.size()) == _V);
	}
	bool startsWith(const char * _V) const
	{
		return (substr(0,strlen(_V)) == _V);
	}
	bool startsWith(const std::vector<std::string> strings) const
	{
		for (const auto & str : strings)
		{
			if (startsWith(str))
				return true;
		}
		return false;
	}
	bool startsWith(const char * str, const char *& p) const
	{
		if (startsWith(str))
		{
			p = c_str() + strlen(str);
			return true;
		}
		else
		{
			p = 0;
			return false;
		}
	};

	void truncateFrom(const char * _V)
	{
		size_t p = find(_V);
		if (p != std::string::npos)
			*this = substr(0,p);
	}
	stringEx & replace (const std::string & str1,const std::string & str2)
	{
		auto i = find(str1);
		if (i == std::string::npos) return *this;
		std::string::replace(i,str1.size(),str2);
		return *this;
	}
	stringEx & replaceAll (const char c1,const char c2)
	{
		size_t _f = 0;
		while ((_f = find_first_of(c1,_f)) != std::string::npos)
			at(_f) = c2;
		return *this;
	}
	stringEx & replaceAll(const std::string & str1, const std::string & str2)
	{
		size_t i;
		while ((i= find(str1)) != npos)
			std::string::replace(i, str1.size(), str2);
		return *this;
	}

	const stringEx operator *(size_t n) const
	{
		stringEx a;
		for (size_t i = 0;i < n;i++)
			a += *this;
		return a;
	}
	const stringEx operator *(int n) const
	{
		stringEx a;
		for (size_t i = 0; i < n; i++)
			a += *this;
		return a;
	}
};

typedef stringEx _s;

inline bool startsWith(const std::string & _S1,const std::string & _S2)
{
	return (_S1.substr(0,_S2.size()) == _S2);
}
inline bool startsWith(const std::string & _S1,const char * _S2)
{
	return (_S1.substr(0,strlen(_S2)) == _S2);
}
inline void toUpper(std::string & _S1)
{
    for(auto& x: _S1)
        x = toupper(x);
}


//	Use fputs so there is no appended \n.  This then reproduces fputs
template <typename... Params>
int printf(const Params& ... params)
{
	return fputs(stringEx(params...).c_str(), stdout);
};



#endif
