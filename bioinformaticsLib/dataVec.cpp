// ***************************************************************************
// DataVec.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// A  vector of values with some standard operators defined
// ***************************************************************************

#include "dataVec.h"
#include "containerEx.h"
#include "parser.h"
#include "rand.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>
#include <mutex>
#include <map>

using namespace std;
//#define DATAVEC_CACHE


#undef dataVec

#ifdef DATAVEC_CACHE

#define CACHE_MINSIZE 500

#ifdef _WIN32
#define THREAD_CACHE  cache[this_thread::get_id()]
map<thread::id,map<size_t,vector<vector<VEC_DATA_TYPE> > > > cache;
#else
thread_local map<size_t,vector<vector<VEC_DATA_TYPE> > > cache;
#define THREAD_CACHE  cache
#endif

dataVec::dataVec(size_t s)
{
	if (s == 0)
		return;
	vector<vector<VEC_DATA_TYPE> > & c = THREAD_CACHE[s];

	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();
	}
	else
		vector<VEC_DATA_TYPE>::resize(s);
};

void dataVec::resize(size_t s)
{
	if (s <  size())
	{
		vector<VEC_DATA_TYPE>::resize(s);
	}
	else if (s >  size())
	{
		vector<vector<VEC_DATA_TYPE> > & c = THREAD_CACHE[s];
		if (c.size())
		{
			swap(c.back());
			for (size_t i = 0;i < c.back().size();i++)
				at(i) = c.back().at(i);
			c.pop_back();
		}
		else
			vector<VEC_DATA_TYPE>::resize(s);
	}
}

dataVec::dataVec(size_t s,VEC_DATA_TYPE v)
{
	if (s == 0)
		return;

	vector<vector<VEC_DATA_TYPE> > & c = THREAD_CACHE[s];
	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();
	}

	assign(s,v);
};

dataVec::dataVec(const std::vector<VEC_DATA_TYPE> & a)
{
	size_t s = a.size();
	vector<vector<VEC_DATA_TYPE> > & c = THREAD_CACHE[s];
	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();
		std::vector<VEC_DATA_TYPE>::operator=(a);
	}
	else
		vector<VEC_DATA_TYPE>::operator=(a);
};
dataVec::dataVec(const dataVec & a)
{
	size_t s = a.size();
	vector<vector<VEC_DATA_TYPE> > & c = THREAD_CACHE[s];
	if (c.size() && (s >= CACHE_MINSIZE))
	{
		swap(c.back());
		c.pop_back();

		std::vector<VEC_DATA_TYPE>::operator=(a);

	}
	else
		vector<VEC_DATA_TYPE>::operator=(a);
};


dataVec::~dataVec()
{
	size_t s = size();
	vector<vector<VEC_DATA_TYPE> > & c = THREAD_CACHE[s];
	if (s < CACHE_MINSIZE)
		return;
	c.emplace_back(move(*this));
}

void dataVec::clearCache()
{
#ifdef _WIN32
	cache[this_thread::get_id()].clear();
#endif
};


#else
//dataVec::dataVec(size_t s): std::vector<VEC_DATA_TYPE>(s){};
//dataVec<>::dataVec<>(size_t s, VEC_DATA_TYPE v):std::vector<VEC_DATA_TYPE>(s,v) {};
//dataVec::dataVec(const vector<VEC_DATA_TYPE> & a): std::vector<VEC_DATA_TYPE>(a){};
//dataVec::dataVec(const dataVec & a):vector<VEC_DATA_TYPE>(a) {};
//dataVec<>::~dataVec<>(){};


template<> void dataVec<>::resize(size_t s) {
	vector<VEC_DATA_TYPE>::resize(s);
};

template<> void dataVec<>::clearCache() {};

#endif


//	Cholesky_Decomposition returns the Cholesky Decomposition Matrix. 
template<> dataVec<> dataVec<>::diagchol() const
{

	dataArray L(size());

	//	Initialize and populate matrix L which will be the lower Cholesky
	for (size_t i = 0; i < size(); i++)
		for (size_t j = 0; j < size(); j++)
		{
			double temp = 0;// , temp2 = 0;
			if (i > j)
			{
//				if (j > 0)
//			{
//					for (k = 1; k < j + 1; k++)
//						temp2 += (L[i][k - 1] * L[j][k - 1]);
//				}
//				L[i][j] = (p[i][j] - temp2) / L[j][j];
			}
			else if (i == j)
			{
				for (size_t k = 0; k < i; k++)
					temp += (L[i][k] * L[i][k]);
				L[i][j] = sqrt(at(i) - temp);
			}
//			else
//				L[i][j] = 0;
		}
	dataVec M(size());
	for (size_t i = 0; i < size(); i++)
		M[i] = L[i][i];	

	return M;
};

template<> dataVec<> dataVec<>::addNoise(double mean) const
{
	size_t s(std::vector<VEC_DATA_TYPE>::size());
	return *this + randn(s)*mean;
}


template<> dataVec<> & dataVec<>::clipMax(double maxVal)
{
	for (auto & i: *this)
		if (i > maxVal)
			i = maxVal;
	return *this;
}




