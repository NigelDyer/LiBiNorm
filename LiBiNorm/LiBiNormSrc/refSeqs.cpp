#include "refSeqs.h"



using namespace std;

bool refSeqs::get(const std::string& ebwtFileBase, std::map<std::string, geneData> & genes)
{
	try
	{
		getInternal<uint32_t>(ebwtFileBase, genes);
		return true;
	}
	catch (...) {};
/*	try
	{
		getInternal<uint64_t>(ebwtFileBase, genes);
		return true;
	}
	catch (...) {};
*/
	return false;
}

