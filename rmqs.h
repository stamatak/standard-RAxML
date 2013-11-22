#ifndef _rmqs_h_
#define _rmqs_h_

#define MEM_COUNT

#include "rmq.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>



static DT *a = NULL;

// size of array a
static DTidx n;

// table M for the out-of-block queries (contains indices of block-minima)
static DTsucc** M = NULL;

// because M just stores offsets (rel. to start of block), this method
// re-calculates the true index:
static inline DTidx m(DTidx k, DTidx block);

// depth of table M:
static DTidx M_depth;

// table M' for superblock-queries (contains indices of block-minima)
static DTidx** Mprime = NULL;

// depth of table M':
static DTidx Mprime_depth;

// type of blocks
static DTsucc2 *type = NULL;

// precomputed in-block queries
static DTsucc** Prec = NULL;

// microblock size
static DTidx s;

// block size
static DTidx sprime;

// superblock size
static DTidx sprimeprime;

// number of blocks (always n/sprime)
static DTidx nb;

// number of superblocks (always n/sprimeprime)
static DTidx nsb;

// number of microblocks (always n/s)
static DTidx nmb;

// return microblock-number of entry i:
static inline DTidx microblock(DTidx i);

// return block-number of entry i:
static inline DTidx block(DTidx i);

// return superblock-number of entry i:
static inline DTidx superblock(DTidx i);

// precomputed Catalan triangle (17 is enough for 64bit computing):
static const DTidx Catalan[17][17];

// minus infinity (change for 64bit version)
static const DT minus_infinity;

// stuff for clearing the least significant x bits (change for 64-bit computing)
static const DTsucc HighestBitsSet[8];
static DTsucc clearbits(DTsucc, DTidx);

// Least Significant Bits for 8-bit-numbers:
static const char LSBTable256[256];

// return least signigicant bit in constant time (change for 64bit version)
static DTidx lsb(DTsucc);

// the following stuff is for fast base 2 logarithms:
// (currently only implemented for 32 bit numbers)
static const char LogTable256[256];
static DTidx log2fast(DTidx);

// set if input array a is very small, so that naive scanning is better:
static bool ARRAY_VERY_SMALL;

#endif
