/* bits - handle operations on arrays of bits. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef BITS_H
#define BITS_H

#include "localmem.h"

typedef unsigned char Bits;

#define bitToByteSize(bitSize) ((bitSize+7)/8)
/* Convert number of bits to number of bytes needed to store bits. */

Bits *bitAlloc(int bitCount);
/* Allocate bits. */

Bits *bitClone(Bits* orig, int bitCount);
/* Clone bits. */

void bitFree(Bits **pB);
/* Free bits. */

Bits *lmBitAlloc(struct lm *lm,int bitCount);
// Allocate bits.  Must supply local memory.

void bitSetOne(Bits *b, int bitIx);
/* Set a single bit. */

void bitClearOne(Bits *b, int bitIx);
/* Clear a single bit. */

void bitSetRange(Bits *b, int startIx, int bitCount);
/* Set a range of bits. */

boolean bitReadOne(Bits *b, int bitIx);
/* Read a single bit. */

int bitCountRange(Bits *b, int startIx, int bitCount);
/* Count number of bits set in range. */

int bitFindSet(Bits *b, int startIx, int bitCount);
/* Find the index of the the next set bit. */

int bitFindClear(Bits *b, int startIx, int bitCount);
/* Find the index of the the next clear bit. */

extern int bitsInByte[256];
/* Lookup table for how many bits are set in a byte. */

void bitsInByteInit();
/* Initialize bitsInByte array. */

#endif /* BITS_H */

