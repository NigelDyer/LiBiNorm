// ***************************************************************************
// DataVec.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 March 2018
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

inline VEC_DATA_TYPE sqr(VEC_DATA_TYPE a) { return a * a; };

template <typename _Ty= VEC_DATA_TYPE, class _Alloc = std::allocator<_Ty> >
class dataVec: public std::vector<_Ty>
{
	typedef std::vector<_Ty> _vec;
public:

	dataVec(const std::vector<_Ty> & a) :_vec(a) {};
	dataVec(size_t s = 0, _Ty v = 0) : _vec(s,v){};
//	dataVec(iterator start,iterator end) :std::vector<_Ty>(start, end) {};

	dataVec(std::initializer_list<_Ty> a) : _vec(a) {};

//	dataVec(const std::vector<VEC_DATA_TYPE> & a);
//	dataVec(const dataVec & a);

/*	dataVec & operator = (const dataVec & a)
	{
		std::vector<VEC_DATA_TYPE>::operator=(a);
		return This;
	}
	*/
	void resize(size_t s = 0);
	void resize(size_t s,const _Ty value);
	size_t size() const { return _vec::size(); };
	_Ty & at(size_t i){ return _vec::at(i); };
	const _Ty & at(size_t i) const { return _vec::at(i); };

#ifdef RVALUE_CONSTRUCTORS
	dataVec(const dataVec & a) : _vec(a)
	{
	};

	dataVec(dataVec && a) : _vec(move(a))
	{
	};

	dataVec & operator = (dataVec && a)
	{
		_vec::operator = (move(a));
		return This;
	}
	dataVec & operator = (const dataVec & a)
	{
		_vec::operator = (a);
		return This;
	}
	dataVec & operator = (_Ty data)
	{
		for (size_t i = 0; i < size(); i++)
			at(i) = data;
//		m_error = 0;
		return (*this);
	}
#endif

	static void clearCache() {};

	~dataVec() {};

	dataVec operator > (_Ty a) const
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = _vec::at(i)>a?1:0;
		return retVal;
	}
	dataVec operator < (const dataVec & a) const
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = _vec::at(i)<a.at(i)?1:0;
		return retVal;
	}
	dataVec operator < (dataVec && a) const
	{
		size_t s = _vec::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) = _vec::at(i)<a.at(i)?1:0;
		return a;
	}
	dataVec operator < (_Ty a) const
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s; i++)
			retVal.at(i) = _vec::at(i)<a ? 1 : 0;
		return retVal;
	}
	dataVec operator * (const dataVec & a) const
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = _vec::at(i) * a.at(i);
		return retVal;
	}
	dataVec operator * (dataVec && a) const
	{
		size_t s = _vec::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) *= _vec::at(i);
		return a;
	}
	dataVec operator / (const dataVec & a) const 
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = _vec::at(i) / a.at(i);
		return retVal;
	}
	dataVec operator / (dataVec && a) const 
	{
		size_t s = _vec::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) = _vec::at(i) / a.at(i);
		return a;
	}
	dataVec operator + (const dataVec & a) const
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = _vec::at(i) + a.at(i);
		return retVal;
	}
	dataVec operator + (dataVec && a) const
	{
		size_t s = _vec::size();
		for (size_t i = 0; i < s;i++)
			a.at(i) += _vec::at(i);
		return a;
	}

	dataVec operator - (const dataVec & a) const 
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = _vec::at(i) - a.at(i);
		return retVal;
	}
	dataVec operator - () const
	{
		size_t s = _vec::size();
		dataVec retVal(s);
		for (size_t i = 0; i < s;i++)
			retVal.at(i) = -_vec::at(i);
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
			return _vec(this->begin()+a, this->end());
		else
			return _vec(this->begin() + a, this->begin() + a + b);
	}

	dataVec & append(const dataVec a)
	{
		insert(_vec::end(), a.begin(), a.end());
		return This;
	}

	template<class _Ty2>
	dataVec & append(const std::vector< _Ty2> a,size_t size)
	{
		size_t addition = std::min(size, a.size());
		size_t len = _vec::size();
		resize(len + addition);
		for (size_t i = 0; i < addition; i++)
			_vec::at(len + i) = a.at(i);
		return This;
	}
	dataVec & append(size_t N, _Ty value)
	{
		insert(_vec::end(),N,value);
		return This;
	}

	dataVec diagchol() const;

	dataVec & clipMax(double maxVal);

	void log()
	{
		size_t count = size();
		for (size_t i = 0; i < count; i++)
			at(i) = ::log(at(i));
	}
	void smooth(int N = 1)
	{
		size_t count = size();
		for (int j = 0; j < N; j++)
		{
			for (size_t i = 1; i < (count - 1); i++)
				at(i) = (at(i-1) + 2 * at(i) + at(i+1))/4;
		}
	}


	_Ty sum()
	{
		_Ty tot = 0;
		for (auto & i : This)
			tot += i;
		return tot;
	}

	void normalise()
	{
		This *= (size()/sum()) ;
	}

	_Ty average()
	{
		return sum()/size();
	}
	_Ty std_dev()
	{
		size_t count = size();
		_Ty ave = average();

		_Ty sigf = 0;
		for (size_t i = 0; i < count; i++)
		{
			sigf += sqr(at(i) - ave);
		}

		sigf = sqrt(sigf / (count - 1));
		return sigf;

	}
	_Ty info_content()
	{
		_Ty info = 0;
		size_t count = size();
		for (size_t i = 0; i < count; i++)
		{
			if (at(i) > 0)
				info += at(i) * log2(at(i));
		}
		return info;
	}


	dataVec expand (const std::vector<int> & i) const
	{
		size_t s(i.size());
		dataVec retVal(s);
		for (size_t j = 0;j < s;j++)
			retVal.at(j) = _vec::at(i.at(j));
		return retVal;
	}
	dataVec addNoise(double mean) const;
};

//	The dataVec template class has a default type of VEC_DATA_TYPE so dataVec<> is all that 
//	is needed to use it with the default.   This #define removes the need for 
//	the <> so we can forget that it is a template class in normal usage
#define dataVec dataVec<>


//	A two dimensiona array
class dataArray : public std::vector<dataVec>
{
public:
	//	For createing a square array of size s
	dataArray(size_t s): std::vector<dataVec>(s,dataVec(s)){};
	//	For creating an s*t rectangular array
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

//	Lots of standard operators
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
inline dataVec operator / (VEC_DATA_TYPE a, dataVec && b)
{
	for (size_t i = 0; i < b.size(); i++)
		b[i] = a / b[i];
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
inline dataVec operator + (VEC_DATA_TYPE b,const dataVec & a)
{
	size_t s(a.size());
	dataVec retVal(s);
	for (size_t i = 0; i < s; i++)
		retVal.at(i) = a.at(i) + b;
	return retVal;
}
inline dataVec operator + (VEC_DATA_TYPE b, dataVec && a)
{
	for (VEC_DATA_TYPE & i : a)
		i += b;
	return a;
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
/*
//	This is required so that a dataVec can be printed using the printEx classes/methods by casting it to
//	the underlying vector
inline bool printVal(outputDataFile * f,const dataVec & value)
{
	printVal(f,(const std::vector<VEC_DATA_TYPE> &) value);
	return true;
};
*/


//	An ordered list that can be used for calculating medians 
//	Add the values using 'add' then 'median' finds the middle 
//	value of the ordered set
class orderedVec : public std::multiset<VEC_DATA_TYPE>
{
public:
	VEC_DATA_TYPE median() const
	{
		if (size() == 0) return 0;
		auto i = begin();
		std::advance(i, size() / 2);

//		size_t n = 0;
//		for (; n < size() / 2; i++, n++) {};
		return *i;
	}
	void add(VEC_DATA_TYPE v) { emplace(v); };
};







#endif
