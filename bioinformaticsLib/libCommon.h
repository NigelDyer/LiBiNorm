// ***************************************************************************
// libCommon.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some common library functions
// ***************************************************************************

#ifndef LIBCOMMON_H
#define LIBCOMMON_H


#ifdef _WIN32
#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG
#else
#define _ASSERT(A) 
#define _ASSERT_EXPR(A,B)
#endif

#include <ctime>
#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include "limits.h"

#if defined(WIN32) || defined(_WIN32) 
#ifndef strdup
#define strdup(A) _strdup(A)
#endif
#define strncasecmp(A,B,C) _strnicmp(A,B,C)
#endif 

#define This (*this)

inline double log2(double _V) {
	static double ln2 = log(2);
	return log(_V)/ln2;

}

#ifdef _DEBUG
#define _DBG(A) A
#else
#define _DBG(A)
#endif


//Allows optional progressMEssages to be output
extern bool verbose;
extern bool debugPrint;

inline void progMessage()
{
	std::cerr << std::endl;
}

#define MAX_DOUBLE std::numeric_limits<double>::max()
#define MIN_DOUBLE std::numeric_limits<double>::min()

template <typename First, typename... Rest>
inline void progMessage(const First & first, const Rest& ... rest)
{
	std::cerr << first;
	progMessage(rest...);
}

template <typename First, typename... Rest>
inline void optMessage(const First & first, const Rest& ... rest)
{
	if (!verbose)
		return;
	std::cerr << first;
	progMessage(rest...);
}

template <typename First, typename... Rest>
inline void debugMessage(const First & first, const Rest& ... rest)
{
	if (!debugPrint)
		return;
	std::cerr << first;
	progMessage(rest...);
}


inline void exitSuccess()
{
	exit(EXIT_SUCCESS);
}

inline void exitFail()
{
#ifdef PAUSE_ON_EXIT_FAILURE
	std::string test;
	std::cin >> test;
#endif
	exit(EXIT_FAILURE);
}
template <typename... Params>
inline void exitFail(const Params& ... params)
{
	progMessage(params...);
	exitFail();
};


//	For timings

inline clock_t & startTime()
{
	static clock_t C;
	return C;
}


inline void initClock()
{
	startTime() = clock();
}

template <typename... Rest>
inline void elapsedTime(const Rest& ... rest)
{
	clock_t now = clock();
	double elapsed_secs = double(now - startTime()) / CLOCKS_PER_SEC;
	optMessage(rest...,": Elapsed time ",elapsed_secs," seconds");
}

class tempDirectory
{
	static std::string & theDirectory();
	~tempDirectory();
public:
	static std::string get(std::string suggestion);
};



#endif
