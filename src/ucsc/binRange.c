/* binRange Stuff to handle binning - which helps us restrict 
 * our attention to the parts of database that contain info
 * about a particular window on a chromosome. This scheme
 * will work without modification for chromosome sizes up
 * to half a gigaBase.  The finest sized bin is 128k (1<<17).
 * The next coarsest is 8x as big (1<<13).  There's a hierarchy
 * of bins with the chromosome itself being the final bin.
 * Features are put in the finest bin they'll fit in. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "binRange.h"


/* add one new level to get coverage past chrom sizes of 512 Mb
 *	effective limit is now the size of an integer since chrom start
 *	and end coordinates are always being used in int's == 2Gb-1 */
static int binOffsetsExtended[] =
	{4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

static int binOffsets[] = {512+64+8+1, 64+8+1, 8+1, 1, 0};
#define _binFirstShift 17	/* How much to shift to get to finest bin. */
#define _binNextShift 3		/* How much to shift to get to next larger bin. */


static int binFromRangeStandard(int start, int end)
/* Given start,end in chromosome coordinates assign it
 * a bin.   There's a bin for each 128k segment, for each
 * 1M segment, for each 8M segment, for each 64M segment,
 * and for each chromosome (which is assumed to be less than
 * 512M.)  A range goes into the smallest bin it will fit in. */
{
int startBin = start, endBin = end-1, i;
startBin >>= _binFirstShift;
endBin >>= _binFirstShift;
for (i=0; i<ArraySize(binOffsets); ++i)
    {
    if (startBin == endBin)
        return binOffsets[i] + startBin;
    startBin >>= _binNextShift;
    endBin >>= _binNextShift;
    }
errAbort("start %d, end %d out of range in findBin (max is 512M)", start, end);
return 0;
}

static int binFromRangeExtended(int start, int end)
/* Given start,end in chromosome coordinates assign it
 * a bin.   There's a bin for each 128k segment, for each
 * 1M segment, for each 8M segment, for each 64M segment,
 * for each 512M segment, and one top level bin for 4Gb.
 *	Note, since start and end are int's, the practical limit
 *	is up to 2Gb-1, and thus, only four result bins on the second
 *	level.
 * A range goes into the smallest bin it will fit in. */
{
int startBin = start, endBin = end-1, i;
startBin >>= _binFirstShift;
endBin >>= _binFirstShift;
for (i=0; i<ArraySize(binOffsetsExtended); ++i)
    {
    if (startBin == endBin)
	return _binOffsetOldToExtended + binOffsetsExtended[i] + startBin;
    startBin >>= _binNextShift;
    endBin >>= _binNextShift;
    }
errAbort("start %d, end %d out of range in findBin (max is 2Gb)", start, end);
return 0;
}

static int binFromRangeBinKeeperExtended(int start, int end)
/* This is just like binFromRangeExtended() above, but it doesn't limit
 * the answers to the range from _binOffsetOldToExtended and up.
 *	It simply uses the whole new bin scheme as if it was the only
 *	one.
 */
{
int startBin = start, endBin = end-1, i;
startBin >>= _binFirstShift;
endBin >>= _binFirstShift;
for (i=0; i<ArraySize(binOffsetsExtended); ++i)
    {
    if (startBin == endBin)
	return binOffsetsExtended[i] + startBin;
    startBin >>= _binNextShift;
    endBin >>= _binNextShift;
    }
errAbort("start %d, end %d out of range in findBin (max is 2Gb)", start, end);
return 0;
}

struct binKeeper *binKeeperNew(int minPos, int maxPos)
/* Create new binKeeper that can cover range. */
{
int binCount;
struct binKeeper *bk;
if (minPos < 0 || maxPos < 0 || minPos > maxPos)
    errAbort("bad range %d,%d in binKeeperNew", minPos, maxPos);

binCount = binFromRangeBinKeeperExtended(maxPos-1, maxPos) + 1;
AllocVar(bk);
bk->minPos = minPos;
bk->maxPos = maxPos;
bk->binCount = binCount;
AllocArray(bk->binLists, binCount);
return bk;
}


void binKeeperAdd(struct binKeeper *bk, int start, int end, void *val)
/* Add item to binKeeper. */ 
{
int bin;
struct binElement *el;
if (start < bk->minPos || end > bk->maxPos || start > end)
    errAbort("(%d %d) out of range (%d %d) in binKeeperAdd", 
    	start, end, bk->minPos, bk->maxPos);
bin = binFromRangeBinKeeperExtended(start, end);
assert(bin < bk->binCount);
AllocVar(el);
el->start = start;
el->end = end;
el->val = val;
slAddHead(&bk->binLists[bin], el);
}

int binElementCmpStart(const void *va, const void *vb)
/* Compare to sort based on start. */
{
const struct binElement *a = *((struct binElement **)va);
const struct binElement *b = *((struct binElement **)vb);
return a->start - b->start;
}

struct binElement *binKeeperFind(struct binKeeper *bk, int start, int end)
/* Return a list of all items in binKeeper that intersect range.
 * Free this list with slFreeList. */
{
struct binElement *list = NULL, *newEl, *el;
int startBin, endBin;
int i,j;

if (start < bk->minPos) start = bk->minPos;
if (end > bk->maxPos) end = bk->maxPos;
if (start >= end) return NULL;
startBin = (start>>_binFirstShift);
endBin = ((end-1)>>_binFirstShift);
for (i=0; i<ArraySize(binOffsetsExtended); ++i)
    {
    int offset = binOffsetsExtended[i];
    for (j=startBin+offset; j<=endBin+offset; ++j)
        {
	for (el=bk->binLists[j]; el != NULL; el = el->next)
	    {
	    if (rangeIntersection(el->start, el->end, start, end) > 0)
	        {
		newEl = CloneVar(el);
		slAddHead(&list, newEl);
		}
	    }
	}
    startBin >>= _binNextShift;
    endBin >>= _binNextShift;
    }
return list;
}

struct binElement *binKeeperFindSorted(struct binKeeper *bk, int start, int end)
/* Like binKeeperFind, but sort list on start coordinates. */
{
struct binElement *list = binKeeperFind(bk, start, end);
slSort(&list, binElementCmpStart);
return list;
}
