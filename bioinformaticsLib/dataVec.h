// ***************************************************************************
// DataVec.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// A  vector of values with some standard operators defined
// ***************************************************************************

#ifndef DATAVEC_H
#define DATAVEC_H

#include <stdlib.h>
#include <vector>
#include <string>
#include <valarray>
#include <algorithm>
#include <set>

#include "libCommon.h"
#include "printEx.h"

#define VEC_DATA_TYPE double
#define RVALUE_CONSTRUCTORS

template <typename _Ty= VEC_DATA_TYPE, class _Alloc = std::allocator<_Ty> >
class dataVec: public std::vector<_Ty>
{
public:

	dataVec(const std::vector<_Ty> & a) :std::vector<_Ty>(a) {};
	dataVec(size_t s = 0, _Ty v = 0) :std::vector<_Ty>(s,v){};
//	dataVec(iterator start,iterator end) :std::vector<_Ty>(start, end) {};

	dataVec(std::initializer_list<_Ty> a) : std::vector<_Ty>(a) {};

//	dataVec(const std::vector<VEC_DATA_TYPE> & a);
//	dataVec(const dataVec & a);

/*	dataVec & operator = (const dataVec & a)
	{
		std::vector<VEC_DATA_TYPE>::operator=(a);
		return This;
	}
	*/
	void resize(size_t s = 0);

#ifdef RVALUE_CONSTRUCTORS
	dataVec(const dataVec & a) : std::vector<_Ty>(a)
	{
	};

	dataVec(dataVec && a) : std::vector<_Ty>(move(a))
	{
	};

	dataVec & operator = (dataVec && a)
	{
		std::vector<_Ty>::operator = (move(a));
		return This;
	}
	dataVec & operator = (const dataVec & a)
	{
		std::vector<_Ty>::operator = (a);
		return This;
	}
#endif

	static void clearCache() {};

	~dataVec() {};

	dataVec operator > (_Ty a) const
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = std::vector<_Ty>::at(i)>a?1:0;
		return retVal;
	}
	dataVec operator < (const dataVec & a) const
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = std::vector<_Ty>::at(i)<a.at(i)?1:0;
		return retVal;
	}
	dataVec operator < (dataVec && a) const
	{
		size_t s = std::vector<_Ty>::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) = std::vector<_Ty>::at(i)<a.at(i)?1:0;
		return a;
	}
	dataVec operator * (const dataVec & a) const
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = std::vector<_Ty>::at(i) * a.at(i);
		return retVal;
	}
	dataVec operator * (dataVec && a) const
	{
		size_t s = std::vector<_Ty>::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) *= std::vector<_Ty>::at(i);
		return a;
	}
	dataVec operator / (const dataVec & a) const 
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = std::vector<_Ty>::at(i) / a.at(i);
		return retVal;
	}
	dataVec operator / (dataVec && a) const 
	{
		size_t s = std::vector<_Ty>::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) = std::vector<_Ty>::at(i) / a.at(i);
		return a;
	}
	dataVec operator + (const dataVec & a) const
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = std::vector<_Ty>::at(i) + a.at(i);
		return retVal;
	}
	dataVec operator + (dataVec && a) const
	{
		size_t s = std::vector<_Ty>::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) += std::vector<_Ty>::at(i);
		return a;
	}

	dataVec operator - (const dataVec & a) const 
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = std::vector<_Ty>::at(i) - a.at(i);
		return retVal;
	}
	dataVec operator - () const
	{
		size_t s = std::vector<_Ty>::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = -std::vector<_Ty>::at(i);
		return retVal;
	}

	bool does_not_contain_null()
	{
		for (_Ty & i : *this)
			if (i == 0)
				return false;
		return true;
	}

	dataVec subset (int a,int b = -1) const
	{
		if (b == -1)
			return std::vector<_Ty>(this->begin()+a, this->end());
		else
			return std::vector<_Ty>(this->begin() + a, this->begin() + a + b);
	}

	dataVec & append(const dataVec a)
	{
		insert(std::vector<_Ty>::end(), a.begin(), a.end());
		return This;
	}

	template<class _Ty2>
	dataVec & append(const std::vector< _Ty2> a,size_t size)
	{
		size_t addition = std::min(size, a.size());
		size_t len = std::vector<_Ty>::size();
		resize(len + addition);
		for (size_t i = 0; i < addition; i++)
			std::vector<_Ty>::at(len + i) = a.at(i);
		return This;
	}
	dataVec & append(size_t N, _Ty value)
	{
		insert(std::vector<_Ty>::end(),N,value);
		return This;
	}

	dataVec diagchol() const;

	dataVec & clipMax(double maxVal);

	dataVec expand (const std::vector<int> & i) const
	{
		size_t s(i.size());
		dataVec retVal(s);
		for (size_t j = 0;j < s;j++)
			retVal.at(j) = std::vector<_Ty>::at(i.at(j));
		return retVal;
	}
	dataVec addNoise(double mean) const;
};

#define dataVec dataVec<>

class dataArray : public std::vector<dataVec>
{
public:
	dataArray(size_t s): std::vector<dataVec>(s,dataVec(s)){};
	dataArray(size_t s, size_t t) : std::vector<dataVec>(s, dataVec(t)) {};
};



inline VEC_DATA_TYPE sum(const dataVec & a)
{
	VEC_DATA_TYPE retVal = 0;
	for (auto & i : a)
		retVal += i;
	return retVal;
}

/*

This is some code for a fast version of the exponent function that 
was found at.   It did not prove to be make much difference so is not currently being used
https://codingforspeed.com/using-faster-exponential-approximation/

inline
double exp1(double x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}
*/
inline dataVec exp(dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
	{
		i = exp(i);
	}
	return a;
}
inline dataVec log(dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
		i = log(i);
	return a;
}
inline dataVec sqrt(dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
		i = sqrt(i);
	return a;
}
inline dataVec operator * (VEC_DATA_TYPE a, const dataVec & b)
{
	size_t s(b.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = b.at(i) * a;
	return retVal;
}

inline dataVec operator * (VEC_DATA_TYPE a, dataVec && b)
{
	for (VEC_DATA_TYPE & i : b)
		i *= a;
	return b;
}

inline dataVec operator < (VEC_DATA_TYPE a, const dataVec & b)
{
	size_t s(b.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a < b.at(i)?1:0;
	return retVal;
}
inline dataVec operator < (VEC_DATA_TYPE a, dataVec && b)
{
	for (VEC_DATA_TYPE & i : b)
		i = a < i?1:0;
	return b;
}

inline dataVec operator - (VEC_DATA_TYPE a, const dataVec & b)
{
	size_t s(b.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a - b.at(i);
	return retVal;

}
inline dataVec operator - (VEC_DATA_TYPE a, dataVec && b)
{
	for (VEC_DATA_TYPE & i : b)
		i = a - i;
	return b;

}
inline dataVec operator + (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i += b;
	return a;
}
inline dataVec operator + (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) + b;
	return retVal;
}

inline dataVec operator += (dataVec & a, const dataVec & b)
{
	_ASSERT_EXPR(a.size() == b.size(), "/= vector sizes do not match");
	for (size_t i = 0; i < a.size(); i++)
		a.at(i) += b.at(i);
	return a;
}

inline dataVec operator - (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i -= b;
	return a;
}
inline dataVec operator - (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) - b;
	return retVal;
}
inline dataVec operator * (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE &  i : a)
		i *= b;
	return a;
}
inline dataVec operator *= (dataVec & a, VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i *= b;
	return a;
}


inline dataVec operator * (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) * b;
	return retVal;
}

inline dataVec operator / (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE &  i : a)
		i /= b;
	return a;
}
inline dataVec operator /= (dataVec & a, VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE & i : a)
		i /= b;
	return a;
}
inline dataVec operator /= (dataVec & a, const dataVec & b)
{
	_ASSERT_EXPR(a.size() == b.size(), "/= vector sizes do not match");
	for (size_t i = 0;i < a.size();i++)
		a.at(i) /= b.at(i);
	return a;
}
inline dataVec operator / (const dataVec & a ,VEC_DATA_TYPE b)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s;i++)
		retVal.at(i) = a.at(i) / b;
	return retVal;
}

inline dataVec operator ^ (const dataVec & a,int b)
{
	dataVec _ret(a);
	for (int p = 0;p < (b - 1);p++)
		for (VEC_DATA_TYPE & i : _ret)
			i *= i;
	return _ret;
}

inline dataVec operator ^ (dataVec && a,int b)
{
	for (int p = 0;p < (b - 1);p++)
		for (VEC_DATA_TYPE & i : a)
			i *= i;
	return a;
}

inline dataVec operator > (dataVec && a,VEC_DATA_TYPE b)
{
	for (VEC_DATA_TYPE &  i : a)
		i = i>b?1:0;
	return a;
}

inline bool printVal(outputDataFile * f,const dataVec & value)
{
	printVal(f,(const std::vector<VEC_DATA_TYPE> &) value);
	return true;
};



/*	An ordered list that can be used for calculating medians */
class orderedVec : public std::multiset<VEC_DATA_TYPE>
{
public:
	VEC_DATA_TYPE median() const
	{
		if (size() == 0) return 0;
		auto i = begin();
		size_t n = 0;
		for (; n < size() / 2; i++, n++) {};
		return *i;
	}
	void add(VEC_DATA_TYPE v) { emplace(v); };
};







#endif
