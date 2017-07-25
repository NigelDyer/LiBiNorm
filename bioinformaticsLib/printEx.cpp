// ***************************************************************************
// printEx.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// For printing text files
// ***************************************************************************

#include "printEx.h"
#include "stringEx.h"

using namespace std;

outputFile::~outputFile()
{
	close();
}

bool outputFile::open(const char * filename)
{
	if (filename)
		fout = fopen(filename,"wb");
	else
		fout = stdout;
	return (fout != 0);
}

//	Pass a zero length string to open and the standard output is used.
bool outputFile::open(const std::string & filename)
{
	if (filename.size())
		return (open(filename.c_str()));
	else
		return (open(0));
}

void outputFile::flush()
{
	fflush(fout);
}
void outputFile::close()
{
	if (fout)
		fclose(fout);
	fout = 0;
}

int outputFile::printf(const char * line,...)
{
	va_list args;
	va_start(args, line);
	return fprintf(fout,line,args);
};



