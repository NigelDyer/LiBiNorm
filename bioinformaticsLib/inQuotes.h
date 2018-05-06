// ***************************************************************************
// inQuotes.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 28 Feb 2018
// ---------------------------------------------------------------------------
// Some functions for parsing and printing text in quotations
// ***************************************************************************

#ifndef INQUOTES_H
#define INQUOTES_H

/********************************************************************

To put any value A within quotation marks use inQuotes(A).  

This calls the template inQuotes method which then creates the inQuotesPrivate class which is recognised and processed by the associated
PrintVal.   Note that inQuotes is also supported by stringEx, and the printVal makes uses of the stringEx method.  

Using a function inQuotes(..) to create the class avoids the need to explicitly state the template class type when creating the class as it is not needed
by the function which determines in from the parameter type and can then use this to create the templated class

*/

template <class _T>
class inQuotesPrivate
{
template <class _Q>
friend inQuotesPrivate<_Q> inQuotes(const _Q & val);

	inQuotesPrivate(const _T & v):value(v){};
public:
	const _T & value;
};

template <class _Q>
inline	inQuotesPrivate<_Q> inQuotes(const _Q & val)
{
	return inQuotesPrivate<_Q>(val);
};


/***************************************************

User specified format using  'fmt', e.g. 

double d = 4.6666;
file.print(fmt("%6f",d));
*/

template <class _Q>
class fmtPrivate
{
//	Oddly, template<> needs to come before friend.  This can only be constructed using fmt.
template <class _R> friend fmtPrivate<_R> fmt(const char * format,const _R & val);
template <class _R> friend fmtPrivate<_R> $(const char * format, const _R & val);

	fmtPrivate(const char * f,_Q v):format(f),value(v){};
public:
	const char * format;
	_Q value;
};

template <class _R>
inline	fmtPrivate<_R> $(const char * format, const _R & val)
{
	return fmtPrivate<_R>(format, val);
};

template <class _R>
inline	fmtPrivate<_R> fmt(const char * format,const _R & val)
{
	return fmtPrivate<_R>(format,val);
};


#endif
