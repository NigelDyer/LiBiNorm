/*
* Copyright 2015, Daehwan Kim <infphilo@gmail.com>
*
* This file was originally part of HISAT 2.
*
* HISAT 2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* HISAT 2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
*/

// ***************************************************************************
// hisat2Lib.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
//	This code was developed from original code in the HISAT 2 library
//	Code taken from gfm.h and hgfm.h


#ifndef HISAT2_LIB
#define HISAT2_LIB


#include <vector>
#include <string>
#include <iostream>
#include <limits>

using namespace std;

/**
* Extended Burrows-Wheeler transform header.  This together with the
* actual data arrays and other text-specific parameters defined in
* class Ebwt constitute the entire Ebwt.
*/
template <typename index_t>
class GFMParams {

public:
	GFMParams() { }

	GFMParams(
		index_t len,
		index_t gbwtLen,
		index_t numNodes,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		index_t eftabLen,
		bool entireReverse)
	{
		init(len, gbwtLen, numNodes, lineRate, offRate, ftabChars, eftabLen, entireReverse);
	}

	void init(
		index_t len,
		index_t gbwtLen,
		index_t numNodes,
		int32_t lineRate,
		int32_t offRate,
		int32_t ftabChars,
		index_t eftabLen,
		bool entireReverse)
	{
		_entireReverse = entireReverse;
		_linearFM = (len + 1 == gbwtLen || gbwtLen == 0);
		_len = len;
		_gbwtLen = (gbwtLen == 0 ? len + 1 : gbwtLen);
		_numNodes = (numNodes == 0 ? len + 1 : numNodes);
		if (_linearFM) {
			_sz = (len + 3) / 4;
			_gbwtSz = _gbwtLen / 4 + 1;
		}
		else {
			_sz = (len + 1) / 2;
			_gbwtSz = _gbwtLen / 2 + 1;
		}
		_lineRate = lineRate;
		_origOffRate = offRate;
		_offRate = offRate;
		_offMask = std::numeric_limits<index_t>::max() << _offRate;
		_ftabChars = ftabChars;
		_eftabLen = eftabLen;
		_eftabSz = _eftabLen * sizeof(index_t);
		_ftabLen = (1 << (_ftabChars * 2)) + 1;
		_ftabSz = _ftabLen * sizeof(index_t);
		_offsLen = (_numNodes + (1 << _offRate) - 1) >> _offRate;
		_offsSz = _offsLen * sizeof(index_t);
		_lineSz = 1 << _lineRate;
		_sideSz = _lineSz * 1 /* lines per side */;
		if (_linearFM) {
			_sideGbwtSz = _sideSz - (sizeof(index_t) * 4);
			_sideGbwtLen = _sideGbwtSz << 2;
		}
		else {
			_sideGbwtSz = _sideSz - (sizeof(index_t) * 6);
			_sideGbwtLen = _sideGbwtSz << 1;
		}
		_numSides = (_gbwtSz + (_sideGbwtSz)-1) / (_sideGbwtSz);
		_numLines = _numSides * 1 /* lines per side */;
		_gbwtTotLen = _numSides * _sideSz;
		_gbwtTotSz = _gbwtTotLen;
		assert(repOk());
	}

#ifndef NDEBUG
	/// Check that this EbwtParams is internally consistent
	bool repOk() const {
		// assert_gt(_len, 0);
		assert_gt(_lineRate, 3);
		assert_geq(_offRate, 0);
		assert_leq(_ftabChars, 16);
		assert_geq(_ftabChars, 1);
		assert_lt(_lineRate, 32);
		assert_lt(_ftabChars, 32);
		assert_eq(0, _gbwtTotSz % _lineSz);
		return true;
	}
#endif
	index_t  _len;
	index_t  _gbwtLen;
	index_t  _sz;
	index_t  _gbwtSz;
	int32_t  _lineRate;
	int32_t  _origOffRate;
	int32_t  _offRate;
	index_t  _offMask;
	int32_t  _ftabChars;
	index_t  _eftabLen;
	index_t  _eftabSz;
	index_t  _ftabLen;
	index_t  _ftabSz;
	index_t  _offsLen;
	index_t  _offsSz;
	index_t  _lineSz;
	index_t  _sideSz;
	index_t  _sideGbwtSz;
	index_t  _sideGbwtLen;
	index_t  _numSides;
	index_t  _numLines;
	index_t  _gbwtTotLen;
	index_t  _gbwtTotSz;
	bool     _entireReverse;
	bool     _linearFM;
	index_t  _numNodes;
};


class GFMFileOpenException : public std::runtime_error {
public:
	GFMFileOpenException(const std::string& msg = "") :
		std::runtime_error(msg) { }
};
/**
* Flags describing type of Ebwt.
*/
enum GFM_FLAGS {
	GFM_ENTIRE_REV = 4 // true -> reverse Ebwt is the whole
					   // concatenated string reversed, rather than
					   // each stretch reversed
};

template <typename index_t>
void readEbwtRefnamesAndLengths(string instr, vector<string>& refnames, vector<index_t>& refLengths) {

	// _in1 must already be open with the get cursor at the
	// beginning and no error flags set.

	ifstream in;
	string gfm_ext = (sizeof(index_t) == 4) ? "ht2" : "ht21";

	// Initialize our primary and secondary input-stream fields
	in.open((instr + ".1." + gfm_ext).c_str(), ios_base::in | ios::binary);
	if (!in.is_open()) {
		throw GFMFileOpenException("Cannot open file " + instr);
	}

	// Read endianness hints from both streams
	bool switchEndian = false;
	uint32_t one = readU32(in, switchEndian); // 1st word of primary stream
	if (one != 1) {
		switchEndian = true;
	}
	readU32(in, switchEndian); // version

							   // Reads header entries one by one from primary stream
	index_t len = readIndex<index_t>(in, switchEndian);
	index_t gbwtLen = readIndex<index_t>(in, switchEndian);
	index_t numNodes = readIndex<index_t>(in, switchEndian);
	int32_t lineRate = readI32(in, switchEndian);
	/*int32_t  linesPerSide =*/ readI32(in, switchEndian);
	int32_t offRate = readI32(in, switchEndian);
	int32_t ftabChars = readI32(in, switchEndian);
	index_t eftabLen = readIndex<index_t>(in, switchEndian);
	// BTL: chunkRate is now deprecated
	int32_t flags = readI32(in, switchEndian);
	bool entireReverse = false;
	if (flags < 0) {
		entireReverse = (((-flags) & GFM_ENTIRE_REV) != 0);
	}

	// Create a new EbwtParams from the entries read from primary stream
	GFMParams<index_t> gh(len, gbwtLen, numNodes, lineRate, offRate, ftabChars, eftabLen, entireReverse);

	index_t nPat = readIndex<index_t>(in, switchEndian); // nPat

	refLengths.resize(nPat);
	if (switchEndian) {
		for (index_t i = 0; i < nPat; i++) {
			refLengths[i] = readIndex<index_t>(in, switchEndian);
		}
	}
	else {
		in.read((char *)refLengths.data(), nPat * sizeof(index_t));
		size_t r = in.gcount();
		if (r != (nPat * sizeof(index_t))) {
			cerr << "Error reading lengths array: " << r << ", " << nPat * sizeof(index_t) << endl;
			throw 1;
		}
	}


	//	in.seekg(nPat * sizeof(index_t), ios_base::cur); // skip plen

	// Skip rstarts
	index_t nFrag = readIndex<index_t>(in, switchEndian);
	in.seekg(nFrag * sizeof(index_t) * 3, ios_base::cur);

	// Skip ebwt
	in.seekg(gh._gbwtTotLen, ios_base::cur);

	// Skip zOff from primary stream
	index_t numZOffs = readIndex<index_t>(in, switchEndian);
	in.seekg(numZOffs * sizeof(index_t), ios_base::cur);

	// Skip fchr
	in.seekg(5 * sizeof(index_t), ios_base::cur);

	// Skip ftab
	in.seekg(gh._ftabLen * sizeof(index_t), ios_base::cur);

	// Skip eftab
	in.seekg(gh._eftabLen * sizeof(index_t), ios_base::cur);

	// Read reference sequence names from primary index file
	while (true) {
		string line;
		if (in.eof())
			break;
		getline(in, line);
		if (line.empty() || (line[0] == '\0'))
			break;
		refnames.push_back(line);
	}

	// Be kind
	in.clear(); in.seekg(0, ios::beg);
}



#endif
