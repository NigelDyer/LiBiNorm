/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

 // ***************************************************************************
 // Modifications (c) 2018 Nigel Dyer
 // School of Life Sciences, University of Warwick
 // ---------------------------------------------------------------------------

#ifndef REF_READ_H_
#define REF_READ_H_

#include <iostream>
#include <cassert>
#include <string>
#include <ctype.h>
#include <fstream>
#include <stdexcept>
#include "alphabet.h"
#include "word_io.h"
#include "endian_swap.h"

using namespace std;
/**
 * Encapsulates a stretch of the reference containing only unambiguous
 * characters.  From an ordered list of RefRecords, one can (almost)
 * deduce the "shape" of the reference sequences (almost because we
 * lose information about stretches of ambiguous characters at the end
 * of reference sequences).
 */
template <typename TIndexOffU> class RefRecord {
public:
	RefRecord() : off(), len(), first() { }
	RefRecord(TIndexOffU _off, TIndexOffU _len, bool _first) :
		off(_off), len(_len), first(_first)
	{ }

	RefRecord(FILE *in, bool swap) {
		assert(in != NULL);
		if(!fread(&off, sizeof(TIndexOffU), 1, in)) {
			cerr << "Error reading RefRecord offset from FILE" << endl;
			throw 1;
		}
		if(swap) off = endianSwapIndex(off);
		if(!fread(&len, sizeof(TIndexOffU), 1, in)) {
			cerr << "Error reading RefRecord offset from FILE" << endl;
			throw 1;
		}
		if(swap) len = endianSwapIndex(len);
		first = fgetc(in) ? true : false;
	}

	void write(std::ostream& out, bool be) {
		writeIndex<TIndexOffU>(out, off, be);
		writeIndex<TIndexOffU>(out, len, be);
		out.put(first ? 1 : 0);
	}

	TIndexOffU off; /// Offset of the first character in the record
	TIndexOffU len; /// Length of the record
	bool   first; /// Whether this record is the first for a reference sequence
};


#endif /*ndef REF_READ_H_*/
