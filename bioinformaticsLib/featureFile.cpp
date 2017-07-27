// ***************************************************************************
// featureFile.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Some functions for processing gtf and gff files
// ***************************************************************************

#include <fstream>
#include "stringEx.h"
#include "featureFile.h"
#include "parser.h"
#include "printEx.h"

using namespace std;

std::map<std::string, unsigned char> entryTags::idMap ;
std::vector<std::string> entryTags::revIdMap ;


string noSpaces(string S)
{
	size_t i= 0;
	while ((i = S.find_first_of(' ',i)) != string::npos)
		S.replace(i,1,"");
	return S;

}


void entryTags::add(const std::string type1,const std::string val)
{
	//	The types are persisted as unsigned char with each type string being mapped to an unsigned char
	//	the idMap is a map from string to unsigned char and the revIdMap is a reverse map
	map<string, unsigned char>::iterator i = idMap.find(type1);
	if (i == idMap.end())
	{
		revIdMap.push_back(type1);
		i = idMap.emplace(type1,(unsigned char)revIdMap.size()-1).first;
	}
	//	i is now a pointer to the position in the map containing the string.
	emplace_back(entryTag(i->second, val));
}
void entryTags::add(const std::string type1, const std::string type2, const std::string val)
{
	//	The types are persisted as unsigned char with each type string being mapped to an unsigned char
	//	the idMap is a map from string to unsigned char and the revIdMap is a reverse map
	map<string, unsigned char>::iterator i1 = idMap.find(type1);
	if (i1 == idMap.end())
	{
		revIdMap.push_back(type1);
		i1 = idMap.emplace(type1, (unsigned char)revIdMap.size() - 1).first;
	}
	map<string, unsigned char>::iterator i2 = idMap.find(type2);
	if (i2 == idMap.end())
	{
		revIdMap.push_back(type2);
		i2 = idMap.emplace(type2, (unsigned char)revIdMap.size() - 1).first;
	}
	//	i is now a pointer to the position in the map containing the string.
	emplace_back(entryTag(i1->second, i2->second, val));
}

//	Creates an entry from a line of text in the feature file.  If id_attribute is specified then
//	only that attribute is persisted (to save space and increase speed), otherwise all of the attributes are persisted
fileEntry::fileEntry(const string & line, parseTagsFunc parseTags, const std::string & id_attribute) :valid(false)
{
	specificEntry specificTags(*this, parseTags,id_attribute);
	parseTsv(line, chromosome, source, type, start, finish, dot, strand, qual, specificTags);
}

fileEntry::fileEntry(const string & chromosome,const string & featType,const featureSection & feature) :
		source("gbFile"),chromosome(chromosome),dot('.'),valid(false),qual(1000)
{

	const string & name = feature["gene"];
	string transcript_id = feature["transcript_id"];
	if (transcript_id.empty())
		transcript_id = name;

	tags.add("gene_id",noSpaces(name));
	tags.add("transcript_id",noSpaces(transcript_id));
	if (featType == "mRNA")
		tags.add("exon_number",stringEx(1));

//	description = stringEx("gene_id ",inQuotes(noSpaces(name)),"; transcript_id ",inQuotes(noSpaces(transcript_id)),";",(featType == "mRNA")?"exon_number 1":"");

	start = feature.start.value;
	finish = feature.finish.value;
	strand = '+';
	type = featType;
}


void fileEntry::output(printTagsFunc printTags,TsvFile & file)
{
	if (valid)
	{
		file.printStart(chromosome,source,type,start,finish,qual,strand,dot);
		printTags(tags,&file);
		file.printEnd();
	}
}

void fileEntry::addOffset(size_t begin, size_t end, int newLen)
{
	if (valid)
	{
		bool largeRegion = ((end - start) > 3000);
		int indel = newLen - (int)(end - begin);
		if (finish > end)
			finish += indel;
		else if (finish > begin)
		{
			if (largeRegion)
			{
				if (finish > (begin + newLen))
					valid = false;
			}
			else
				finish += (int)(indel * ((double)(finish - begin)/(end - begin)));
		}

		if (start > end)
			start += indel;
		else if (start > begin)
		{
			if (largeRegion)
			{
				if (finish > (start + newLen))
					valid = false;
			}
			else
				start += (int)(indel * ((double)(start - begin)/(end - begin)));
		}
	}
}

void featureFile::addEntry (fileEntry && entry)
{
	entryMap[entry.chromosome].emplace(entry.start,move(entry));
}

bool featureFile::printEntries(const string & filename)
{
	TsvFile file;
	if (!file.open(filename))
		return false;
	for (auto i = entryMap.begin(); i != entryMap.end();i++)
	{
		for (auto j = i->second.begin(); j != i->second.end(); j++)
		{
			file.printStart(i->first, _s(j->second.chromosome,":", j->second.start,"-",j->second.finish), j->second.strand, j->second.name ,
				j->second.source, j->second.type,j->second.biotype,j->second.qual,j->second.valid);

			printTags(j->second.tags,&file);
			file.printEnd();
		}
	}
	file.close();
	return true;
}

//The tags in gff and gtf files have different formats.  The following are specialisations of the parser and print functions
//	that deal with the differences
//	Gtf tags are the easiest, they are just [<tagType> <tagValue>;]
bool printGtfTags(const entryTags & tags, outputDataFile * f)
{
	for (const auto & t : tags)
		fprintf(f->fout, "%s \"%s\";", entryTags::revIdMap[t.type1].c_str(), t.val.c_str());

	return true;
};

//	
void parseGtfTags(fileEntry & fe, const char * start, size_t & len, const stringEx & type)
{
	std::string type1, val;
	while (*start)
	{
		parser(start, " ", type1, val, "; ");
//		if (!htSeqCompatible && ((type1 == "transcript_biotype") && (val == "retained_intron")))
//			tags.allowed = false;
//		else
		if (!type)
			fe.tags.add(type1, val);
		else if (type1 == type)
			fe.name = val;
		if (type1 == "transcript_biotype")
			fe.biotype = val;
		if ((type1 == "gene_biotype") && (!fe.biotype))	//exons have transcript and gene biotypes.  Transcript biotype takes precedence
			fe.biotype = val;
	}
};

void parseGffTags(fileEntry & fe, const char * start, size_t & len, const stringEx & type)
{
	string type1, type2, val1, val2;
	while (*start)
	{
		parser(start, "=;", type1, val1);
		if (!type) // || (type1 == type) || (type1 == "gene_biotype") || (type1 == "transcript_biotype"))
			fe.tags.add(type1, val1);
		else if (type1 == type)
			fe.name = val1;

		if (type1 == "gene_biotype")
			fe.biotype = val1;
		else if (type1 == "Dbxref")
		{
			const char * p = val1.c_str();
			while (*p)
			{
				parser(p, ":,", type2, val2);
				if (!type)
					fe.tags.add(type2, type1, val2);
				else if (type2 == type)
					fe.name = val2;
				if (type2 == "GeneID")
					fe.id = val2;
			}
		}
	}
};

bool printGffTags(const entryTags & tags,outputDataFile * f)
{
	bool semiColonNeeded = false;
	char multiSet = -1;
	for (const auto & t : tags)
	{
		if (t.type2 != multiSet)
		{
			if (semiColonNeeded)fprintf(f->fout, ";");
			if (t.type2 != -1)
				fprintf(f->fout, "%s=", entryTags::revIdMap[t.type2].c_str());

		}
		else
		{
			if (semiColonNeeded)fprintf(f->fout, (t.type2 == -1) ? ";" : ",");
		}
		semiColonNeeded = true;
		char sep2((t.type2 == -1) ? '=' : ':');
		fprintf(f->fout, "%s%c%s", entryTags::revIdMap[t.type1].c_str(), sep2, t.val.c_str());
	}

	return true;
};
/*
setEx<string> allowed = { { "protein_coding",
	"TEC",
	"IG_C_gene","IG_C_pseudogene","IG_J_gene","IG_D_gene","IG_D_pseudogene","IG_V_pseudogene","IG_V_gene","IG_LV_gene",//OK
	"TR_V_gene","TR_J_gene","TR_V_pseudogene","TR_D_gene","TR_V_pseudogene","TR_J_pseudogene","TR_C_gene",//OK
	"retained_intron", 
	"sense_intronic","non_stop_decay", //OK
	"pseudogene","polymorphic_pseudogene",	"transcribed_unprocessed_pseudogene","translated_processed_pseudogene","unitary_pseudogene", //OK
	"transcribed_unprocessed_pseudogene","transcribed_processed_pseudogene","nonsense_mediated_decay","unprocessed_pseudogene",//OK
	"processed_pseudogene","translated_unprocessed_pseudogene",//OK
	"3prime_overlapping_ncrna","lincRNA","macro_lncRNA",//OK
	"miRNA","misc_RNA",//OK
	"sense_overlapping","antisense",//OK
	"Mt_tRNA","Mt_rRNA","rRNA",//OK
	"snRNA","snoRNA","sRNA",//OK
	"ribozyme","scaRNA",//OK
	"processed_transcript",//OK
	} };
*/

//	Opens a gff or gtf format feature file.   
bool featureFile::open(const stringEx & filename,const string & id_attribute,const string & feature_type, bool mapFile)
{
	size_t lineCounter(0);
	ifstream file;
	bool gffFile = false;
	file.open(filename);
	if (!file.is_open()) return false;

	TsvFile mappedFile;
	stringEx line,currChromosome;
	_DBG(int count = 0;)
	getline(file,line);
	if (line.startsWith("##gff-version"))
	{
		if (mapFile)
		{
			//	Mapped feature file is sent to stdout
			mappedFile.open();
			if (!mappedFile.is_open())
				exitFail("Failed to open stdout for mapped feature file");
			mappedFile.print(line);
		}
		gffFile = true;
		getline(file, line);
		//	Setup the parse and print methods that are specific to the feature file type
		parseTags = parseGffTags;
		printTags = printGffTags;
	}
	else
	{
		parseTags = parseGtfTags;
		printTags = printGtfTags;
	}

	bool useChromosome = true;

	map<string,string> geneTypes;
	map<string, string> geneNames;

	do {
		if (line[0] == '#')
		{
			if (line.startsWith("#!"))
			{
				headers.push_back(line);
				if (mapFile) mappedFile.print(line);
			}
			else if (line.startsWith("##"))			//gff specific data
			{
				if (line.startsWith("##sequence-region"))
				{
					string params[2];
					int start, finish;
					parser(line, " \n\r", params, start, finish);

					auto i = chromosomeMap.find(finish);
					if (i != chromosomeMap.end())
					{
						useChromosome = true;
						currChromosome = i->second;

						optMessage(params[1]," maps to chromosome ",currChromosome);
						if (mapFile) mappedFile.print(params[0],currChromosome,start,finish);
					}
					else
					{
						useChromosome = false;
						progMessage("No mapping for chromsome ",params[1]);
						currChromosome = params[1];
					}
				}
				else if (mapFile) mappedFile.print(line);
			}
			else if (mapFile) mappedFile.print(line);
		}
		else
		{
			if (useChromosome)
			{
				if (mappedFile.is_open())
				{
					fileEntry gtf(line, parseTags);
					if ((gtf.type != "region") && (gtf.source == "BestRefSeq"))
					{
						gtf.chromosome = currChromosome;
						mappedFile.print(fileEntryForPrinting(gtf, printTags));
					}
				}
				else
				{
					//	If the id_attribute of one or more feature types have been specified then just persist
					//	the information that has been requested
					if ((id_attribute.size()) || (feature_type.size()))
					{
						fileEntry fe(line, parseTags, id_attribute);
						if (fe.type == "gene")
						{
							if (!gffFile || fe.source.startsWith("BestRefSeq"))
							{
								geneTypes[fe.id] = fe.biotype;
								geneNames[fe.id] = fe.name;
							}
						}
						else if (fe.type == feature_type)
						{
							if (gffFile)
							{
								if (fe.source == "BestRefSeq")
								{
									fe.chromosome = currChromosome;
									if (!fe.biotype)
										fe.biotype = geneTypes[fe.id];
									if (!fe.name)
										fe.name = geneNames[fe.id];
									//	Only save it if it has a name
									if (fe.name)
										addEntry(move(fe));
								}
							}
							else
							{
								if (!ignoredTranscriptTypes.contains(fe.biotype))
								{
									if (strncasecmp(fe.chromosome.c_str(), "chr", 3) == 0)
										fe.chromosome = fe.chromosome.substr(3);
									//	Once we have used the transcript type, replace it with the gene type
									fe.biotype = geneTypes[fe.id];
									if (fe.name)
										addEntry(move(fe));
								}
							}
						}
					}
					else
					{
						//	Otherwise persist the whole line.
						addEntry(fileEntry(line, parseTags));
					}
				}
			}
		}
		if ((++lineCounter % 100000) == 0)
			optMessage(lineCounter," Feature file lines processed.");

	} while (getline(file,line) _DBG (&& (count++ < 100000000)));
	optMessage(lineCounter," Feature file lines processed.");
	file.close();
	return true;
}
void featureFile::addGene(const string & chromosome,geneInfo & I)
{
	bool added = false;
	for (featureSection & f : I.features)
	{
		if (f.type == "mRNA")
		{
			added = true;
			addEntry(fileEntry(chromosome,"transcript",f));
			addEntry(fileEntry(chromosome,"exon",f));

			if (f.linkedCDS != -1)
			{
				addEntry(fileEntry(chromosome,"CDS",I.features[f.linkedCDS]));
			}
		}
	}
	if (!added)
	{
		for (featureSection & f : I.features)
		{
			if (f.type == "gene")
			{
				addEntry(fileEntry(chromosome,"gene",f));
			}
			else if (f.type == "CDS")
			{
				addEntry(fileEntry(chromosome,"transcript",f));
				addEntry(fileEntry(chromosome,"CDS",f));
			}
		}
	}
}

void featureFile::addOffset(const string & chromosome,int begin,int end, int newLen)
{
	for (auto & entry : entryMap[chromosome])
		entry.second.addOffset(begin,end,newLen);
}


bool featureFile::save(const string & filename)
{
	TsvFile file;
	if (!file.open(filename))
		return false;

	if (headers.empty())
		file.print("track name=genbank_transcriptome description=genbank_transcriptome");
	else
	{
		for (auto & h : headers)
			file.print(h);
	}

	for (auto & entries : entryMap)
		for (auto & entry : entries.second)
			entry.second.output(printTags,file);

	file.close();
	return true;
}

//	Some additional methods that are picked up by the parser class and used for
//	feature file specific classes
namespace parserInternal
{
	void parseval(const char * start, specificEntry & value, size_t & len)
	{
		value.parseTags(value.fe,start, len, value.type);
	};

	void parseval(const char * start, qualType & value, size_t & len)
	{
		if ((len == 1) && (*start == '.'))
			value.val = -1;
		else
			parseval(start, value.val, len);
	}
}

