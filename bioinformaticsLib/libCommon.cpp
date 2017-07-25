// ***************************************************************************
// libCommon.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some common library functions
// ***************************************************************************

#ifdef _WIN32
#include <direct.h>
#else
//	For rmdir
#include <unistd.h>	
//	for mkdir
#include <sys/stat.h>
//#define mkdir(A) mkdir(A,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#define mkdir(A) mkdir(A,S_IRWXU )
#endif

#include <chrono>
#include "libCommon.h"
#include "stringEx.h"

using namespace std;
//Controls output of optional progress Messages
bool verbose = false;
bool debugPrint = false;


string & tempDirectory::theDirectory()
{
	static string dir;
	return dir;
}
string tempDirectory::get(string suggestion)
{

	//Only need temporary directory for caching files
	// if the data are position ordered
	//	If a count filename is explicitly specified then we can use this as the basis for the temp directory name
	if (theDirectory().size())
	{
		if (mkdir(theDirectory().c_str()) != 0)
		{
			rmdir(theDirectory().c_str());
			if (mkdir(theDirectory().c_str()) != 0)
				exitFail("Unable create temporary directory ", theDirectory(), "\n Try using the -c option for the count files instead");
		}
	}
	else
	{
		//	Otherwise, put one in the temporary directory
		const char * p = getenv("TEMP");
		if (p == 0)
			p = getenv("TMPDIR");

		if (p == 0)
			theDirectory() = "/tmp/";
		else
			theDirectory() = p;

		srand(time(NULL));
		if ((theDirectory()[theDirectory().size() - 1] != '/') &&
			(theDirectory()[theDirectory().size() - 1] != '\\'))
			theDirectory() += "/";


		auto now_us = chrono::time_point_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now());
		unsigned randValue = chrono::duration_cast<chrono::microseconds>(now_us.time_since_epoch()).count();

		string tempDirRoot = theDirectory();
		theDirectory() += _s("LiBiNorm_temp_", randValue);

		optMessage("temp Directory = ", theDirectory());

		if (mkdir(theDirectory().c_str()) != 0)
		{
			randValue++;

			theDirectory() = _s(tempDirRoot, "LiBiNorm_temp_", randValue);

			optMessage("second attempt at temp Directory = ", theDirectory());

			if (mkdir(theDirectory().c_str()) != 0)
				exitFail("Unable create temporary directory ", theDirectory(), "\n Try using the -c option for the count files instead");
		}
	}
	theDirectory() += "/";
	return theDirectory();
}
tempDirectory::~tempDirectory()
{
	if (theDirectory().size())
	{
		rmdir(theDirectory().c_str());
	}
}


