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

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <fcntl.h>
#include <sys/stat.h>
#include <utility>
#include  <memory.h>
#include "ref_read.h"


/**
 * Concrete reference representation that bulk-loads the reference from
 * the bit-pair-compacted binary file and stores it in memory also in
 * bit-pair-compacted format.  The user may request reference
 * characters either on a per-character bases or by "stretch" using
 * getBase(...) and getStretch(...) respectively.
 *
 * Most of the complexity in this class is due to the fact that we want
 * to represent references with ambiguous (non-A/C/G/T) characters but
 * we don't want to use more than two bits per base.  This means we
 * need a way to encode the ambiguous stretches of the reference in a
 * way that is external to the bitpair sequence.  To accomplish this,
 * we use the RefRecords vector, which is stored in the .3.ebwt index
 * file.  The bitpairs themselves are stored in the .4.ebwt index file.
 *
 * Once it has been loaded, a BitPairReference is read-only, and is
 * safe for many threads to access at once.
 */
template <class TIndexOffU> class BitPairReference {

public:
	/**
	 * Load from .3.ebwt/.4.ebwt Bowtie index files.
	 */
	 /*	BitPairReference(
			 const string& in,
			 bool sanity = false,
			 bool verbose = false,
			 bool startVerbose = false);
	~BitPairReference();
	*/

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	int getStretchNaive(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count) const;

	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	int getStretch(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count
		ASSERT_ONLY(, vector<uint32_t>& destU32_2)) const;
	*/
	/**
	 * Return the number of reference sequences.
	 */
	TIndexOffU numRefs() const {
		return nrefs_;
	}
protected:

	uint32_t byteToU32_[256];

	vector<RefRecord<TIndexOffU>> recs_;       /// records describing unambiguous stretches
	// following two lists are purely for the binary search in getStretch
	vector<TIndexOffU> cumUnambig_; // # unambig ref chars up to each record
	vector<TIndexOffU> cumRefOff_;  // # ref chars up to each record
	vector<TIndexOffU> refLens_;    /// approx lens of ref seqs (excludes trailing ambig chars)
	vector<TIndexOffU> refOffs_;    /// buf_ begin offsets per ref seq
	vector<TIndexOffU> refRecOffs_; /// record begin/end offsets per ref seq
	uint8_t *buf_;      /// the whole reference as a big bitpacked byte array
	TIndexOffU bufSz_;    /// size of buf_
	TIndexOffU bufAllocSz_;
	TIndexOffU nrefs_;    /// the number of reference sequences
	bool     loaded_;   /// whether it's loaded
	bool     sanity_;   /// do sanity checking
	bool     verbose_;
	ASSERT_ONLY(vector<uint32_t> tmp_destU32_);

public:

	BitPairReference(
		const string& in,
		bool sanity,
		bool verbose,
		bool startVerbose) :
		buf_(NULL),
		loaded_(true),
		sanity_(sanity),
		verbose_(verbose)
	{
		string gfm_ext = (sizeof(TIndexOffU) == 4) ? "ht2" : "ht21";
		string s3 = in + ".3." + gfm_ext;
		string s4 = in + ".4." + gfm_ext;

		FILE *f3, *f4;
		if ((f3 = fopen(s3.c_str(), "rb")) == NULL) {
			cerr << "Could not open reference-string index file " << s3 << " for reading." << endl;
			cerr << "This is most likely because your index was built with an older version" << endl
				<< "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
				<< "index (or download one from the Bowtie website) and try again." << endl;
			loaded_ = false;
			return;
		}
		if ((f4 = fopen(s4.c_str(), "rb")) == NULL) {
			cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
			loaded_ = false;
			return;
		}

		// Read endianness sentinel, set 'swap'
		uint32_t one;
		bool swap = false;
		one = readIndex<int32_t>(f3, swap);
		if (one != 1) {
			assert_eq(0x1000000, one);
			swap = true; // have to endian swap U32s
		}

		// Read # records
		TIndexOffU sz;
		sz = readIndex<TIndexOffU>(f3, swap);
		if (sz == 0) {
			cerr << "Error: number of reference records is 0 in " << s3.c_str() << endl;
			throw 1;
		}

		// Read records
		nrefs_ = 0;

		// Cumulative count of all unambiguous characters on a per-
		// stretch 8-bit alignment (i.e. count of bytes we need to
		// allocate in buf_)
		TIndexOffU cumsz = 0;
		TIndexOffU cumlen = 0;
		// For each unambiguous stretch...
		for (TIndexOffU i = 0; i < sz; i++) {
			recs_.push_back(RefRecord<TIndexOffU>(f3, swap));
			if (recs_.back().first) {
				// This is the first record for this reference sequence (and the
				// last record for the one before)
				refRecOffs_.push_back((TIndexOffU)recs_.size() - 1);
				// refOffs_ links each reference sequence with the total number of
				// unambiguous characters preceding it in the pasted reference
				refOffs_.push_back(cumsz);
				if (nrefs_ > 0) {
					// refLens_ links each reference sequence with the total number
					// of ambiguous and unambiguous characters in it.
					refLens_.push_back(cumlen);
				}
				cumlen = 0;
				nrefs_++;
			}
			else if (i == 0) {
				cerr << "First record in reference index file was not marked as "
					<< "'first'" << endl;
				throw 1;
			}
			cumUnambig_.push_back(cumsz);
			cumRefOff_.push_back(cumlen);
			cumsz += recs_.back().len;
			cumlen += recs_.back().off;
			cumlen += recs_.back().len;
		}
		if (verbose_ || startVerbose) {
			cerr << "Read " << nrefs_ << " reference strings from "
				<< sz << " records: ";
		}
		// Store a cap entry for the end of the last reference seq
		refRecOffs_.push_back((TIndexOffU)recs_.size());
		refOffs_.push_back(cumsz);
		refLens_.push_back(cumlen);
		cumUnambig_.push_back(cumsz);
		cumRefOff_.push_back(cumlen);
		bufSz_ = cumsz;
		assert_eq(nrefs_, refLens_.size());
		assert_eq(sz, recs_.size());
		if (f3 != NULL) fclose(f3); // done with .3.gfm_ext file
									// Round cumsz up to nearest byte boundary
		if ((cumsz & 3) != 0) {
			cumsz += (4 - (cumsz & 3));
		}
		bufAllocSz_ = cumsz >> 2;
		assert_eq(0, cumsz & 3); // should be rounded up to nearest 4

								 // Allocate a buffer to hold the reference string
		try {
			buf_ = new uint8_t[cumsz >> 2];
			if (buf_ == NULL) throw std::bad_alloc();
		}
		catch (std::bad_alloc) {
			cerr << "Error: Ran out of memory allocating space for the bitpacked reference.  Please" << endl
				<< "re-run on a computer with more memory." << endl;
			throw 1;
		}

		// Open the bitpair-encoded reference file
		f4 = fopen(s4.c_str(), "rb");
		if (f4 == NULL) {
			cerr << "Could not open reference-string index file " << s4.c_str() << " for reading." << endl;
			cerr << "This is most likely because your index was built with an older version" << endl
				<< "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
				<< "index (or download one from the Bowtie website) and try again." << endl;
			loaded_ = false;
			return;
		}
		// Read the whole thing in
		size_t ret = fread(buf_, 1, cumsz >> 2, f4);
		// Didn't read all of it?
		if (ret != (cumsz >> 2)) {
			cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4.c_str() << endl;
			throw 1;
		}
		// Make sure there's no more
		char c;
		ret = fread(&c, 1, 1, f4);
		assert_eq(0, ret); // should have failed
		fclose(f4);

		// Populate byteToU32_
		bool big = currentlyBigEndian();
		for (int i = 0; i < 256; i++) {
			uint32_t word = 0;
			if (big) {
				word |= ((i >> 0) & 3) << 24;
				word |= ((i >> 2) & 3) << 16;
				word |= ((i >> 4) & 3) << 8;
				word |= ((i >> 6) & 3) << 0;
			}
			else {
				word |= ((i >> 0) & 3) << 0;
				word |= ((i >> 2) & 3) << 8;
				word |= ((i >> 4) & 3) << 16;
				word |= ((i >> 6) & 3) << 24;
			}
			byteToU32_[i] = word;
		}
	}


	~BitPairReference() {
		if (buf_ != NULL) delete[] buf_;
	}

	/**
	* Load a stretch of the reference string into memory at 'dest'.
	*
	* This implementation scans linearly through the records for the
	* unambiguous stretches of the target reference sequence.  When
	* there are many records, binary search would be more appropriate.
	*/
	int getStretchNaive(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count) const
	{
		uint8_t *dest = (uint8_t*)destU32;
		uint64_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint64_t recf = refRecOffs_[tidx + 1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint64_t cur = 0;
		uint64_t bufOff = refOffs_[tidx];
		uint64_t off = 0;
		// For all records pertaining to the target reference sequence...
		for (uint64_t i = reci; i < recf; i++) {
			assert_geq(toff, off);
			off += recs_[i].off;
			for (; toff < off && count > 0; toff++) {
				dest[cur++] = 4;
				count--;
			}
			if (count == 0) break;
			assert_geq(toff, off);
			if (toff < off + recs_[i].len) {
				bufOff += (TIndexOffU)(toff - off); // move bufOff pointer forward
			}
			else {
				bufOff += recs_[i].len;
			}
			off += recs_[i].len;
			for (; toff < off && count > 0; toff++) {
				assert_lt(bufOff, bufSz_);
				const uint64_t bufElt = (bufOff) >> 2;
				const uint64_t shift = (bufOff & 3) << 1;
				dest[cur++] = (buf_[bufElt] >> shift) & 3;
				bufOff++;
				count--;
			}
			if (count == 0) break;
			assert_geq(toff, off);
		} // end for loop over records
		  // In any chars are left after scanning all the records,
		  // they must be ambiguous
		while (count > 0) {
			count--;
			dest[cur++] = 4;
		}
		assert_eq(0, count);
		return 0;
	}

	/**
	* Load a stretch of the reference string into memory at 'dest'.
	*/
	int getStretch(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count
		ASSERT_ONLY(, vector<uint32_t>& destU32_2)) const
	{
		ASSERT_ONLY(size_t origCount = count);
		ASSERT_ONLY(size_t origToff = toff);
		if (count == 0) return 0;
		uint8_t *dest = (uint8_t*)destU32;
#ifndef NDEBUG
		destU32_2.clear();
		uint8_t *dest_2 = NULL;
		int off2;
		if ((rand() % 10) == 0) {
			destU32_2.resize((origCount >> 2) + 2);
			off2 = getStretchNaive(destU32_2.data(), tidx, origToff, origCount);
			dest_2 = ((uint8_t*)destU32_2.data()) + off2;
		}
#endif
		destU32[0] = 0x04040404; // Add Ns, which we might end up using later
		uint64_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint64_t recf = refRecOffs_[tidx + 1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint64_t cur = 4; // keep a cushion of 4 bases at the beginning
		uint64_t bufOff = refOffs_[tidx];
		uint64_t off = 0;
		int64_t offset = 4;
		bool firstStretch = true;
		bool binarySearched = false;
		uint64_t left = reci;
		uint64_t right = recf;
		uint64_t mid = 0;
		// For all records pertaining to the target reference sequence...
		for (uint64_t i = reci; i < recf; i++) {
			uint64_t origBufOff = bufOff;
			assert_geq(toff, off);
			if (firstStretch && recf > reci + 16) {
				// binary search finds smallest i s.t. toff >= cumRefOff_[i]
				while (left < right - 1) {
					mid = left + ((right - left) >> 1);
					if (cumRefOff_[mid] <= toff)
						left = mid;
					else
						right = mid;
				}
				off = cumRefOff_[left];
				bufOff = cumUnambig_[left];
				origBufOff = bufOff;
				i = left;
				assert(cumRefOff_[i + 1] == 0 || cumRefOff_[i + 1] > toff);
				binarySearched = true;
			}
			off += recs_[i].off; // skip Ns at beginning of stretch
			assert_gt(count, 0);
			if (toff < off) {
				size_t cpycnt = min((size_t)(off - toff), count);
				memset(&dest[cur], 4, cpycnt);
				count -= cpycnt;
				toff += cpycnt;
				cur += cpycnt;
				if (count == 0) break;
			}
			assert_geq(toff, off);
			if (toff < off + recs_[i].len) {
				bufOff += toff - off; // move bufOff pointer forward
			}
			else {
				bufOff += recs_[i].len;
			}
			off += recs_[i].len;
			assert(off == cumRefOff_[i + 1] || cumRefOff_[i + 1] == 0);
			assert(!binarySearched || toff < off);
			if (toff < off) {
				if (firstStretch) {
					if (toff + 8 < off && count > 8) {
						// We already added some Ns, so we have to do
						// a fixup at the beginning of the buffer so
						// that we can start clobbering at cur >> 2
						if (cur & 3) {
							offset -= (cur & 3);
						}
						uint64_t curU32 = cur >> 2;
						// Do the initial few bases
						if (bufOff & 3) {
							const uint64_t bufElt = (bufOff) >> 2;
							const int64_t low2 = bufOff & 3;
							// Lots of cache misses on the following line
							destU32[curU32] = byteToU32_[buf_[bufElt]];
							for (int j = 0; j < low2; j++) {
								((char *)(&destU32[curU32]))[j] = 4;
							}
							curU32++;
							offset += low2;
							const int64_t chars = 4 - low2;
							count -= chars;
							bufOff += chars;
							toff += chars;
						}
						assert_eq(0, bufOff & 3);
						uint64_t bufOffU32 = bufOff >> 2;
						uint64_t countLim = count >> 2;
						uint64_t offLim = ((off - (toff + 4)) >> 2);
						uint64_t lim = min(countLim, offLim);
						// Do the fast thing for as far as possible
						for (uint64_t j = 0; j < lim; j++) {
							// Lots of cache misses on the following line
							destU32[curU32] = byteToU32_[buf_[bufOffU32++]];
#ifndef NDEBUG
							if (dest_2 != NULL) {
								assert_eq(dest[(curU32 << 2) + 0], dest_2[(curU32 << 2) - offset + 0]);
								assert_eq(dest[(curU32 << 2) + 1], dest_2[(curU32 << 2) - offset + 1]);
								assert_eq(dest[(curU32 << 2) + 2], dest_2[(curU32 << 2) - offset + 2]);
								assert_eq(dest[(curU32 << 2) + 3], dest_2[(curU32 << 2) - offset + 3]);
							}
#endif
							curU32++;
						}
						toff += (lim << 2);
						assert_leq(toff, off);
						assert_leq((lim << 2), count);
						count -= (lim << 2);
						bufOff = bufOffU32 << 2;
						cur = curU32 << 2;
					}
					// Do the slow thing for the rest
					for (; toff < off && count > 0; toff++) {
						assert_lt(bufOff, bufSz_);
						const uint64_t bufElt = (bufOff) >> 2;
						const uint64_t shift = (bufOff & 3) << 1;
						dest[cur++] = (buf_[bufElt] >> shift) & 3;
						bufOff++;
						count--;
					}
					firstStretch = false;
				}
				else {
					// Do the slow thing
					for (; toff < off && count > 0; toff++) {
						assert_lt(bufOff, bufSz_);
						const uint64_t bufElt = (bufOff) >> 2;
						const uint64_t shift = (bufOff & 3) << 1;
						dest[cur++] = (buf_[bufElt] >> shift) & 3;
						bufOff++;
						count--;
					}
				}
			}
			if (count == 0) break;
			assert_eq(recs_[i].len, bufOff - origBufOff);
			assert_geq(toff, off);
		} // end for loop over records
		  // In any chars are left after scanning all the records,
		  // they must be ambiguous
		while (count > 0) {
			count--;
			dest[cur++] = 4;
		}
		assert_eq(0, count);
		return (int)offset;
	}

};

#endif
