/* bigBed - interface to binary file with bed-style values (that is a bunch of
 * possibly overlapping regions. */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "hash.h"
#include "linefile.h"
#include "obscure.h"
#include "dystring.h"
#include "rangeTree.h"
#include "cirTree.h"
#include "bPlusTree.h"
#include "basicBed.h"
#include "asParse.h"
#include "zlibFace.h"
#include "sig.h"
#include "udc.h"
#include "bbiFile.h"
#include "bigBed.h"

struct bbiFile *bigBedFileOpen(char *fileName)
/* Open up big bed file. */
{
return bbiFileOpen(fileName, bigBedSig, "big bed");
}


struct bigBedInterval *bigBedIntervalQuery(struct bbiFile *bbi, char *chrom,
	bits32 start, bits32 end, int maxItems, struct lm *lm)
/* Get data for interval.  Return list allocated out of lm.  Set maxItems to maximum
 * number of items to return, or to 0 for all items. */
{
struct bigBedInterval *el, *list = NULL;
int itemCount = 0;
bbiAttachUnzoomedCir(bbi);
// Find blocks with padded start and end to make sure we include zero-length insertions:
bits32 paddedStart = (start > 0) ? start-1 : start;
bits32 paddedEnd = end+1;
bits32 chromId;
struct fileOffsetSize *blockList = bbiOverlappingBlocks(bbi, bbi->unzoomedCir,
	chrom, paddedStart, paddedEnd, &chromId);
struct fileOffsetSize *block, *beforeGap, *afterGap;
struct udcFile *udc = bbi->udc;
boolean isSwapped = bbi->isSwapped;

/* Set up for uncompression optionally. */
char *uncompressBuf = NULL;
if (bbi->uncompressBufSize > 0)
    uncompressBuf = needLargeMem(bbi->uncompressBufSize);

char *mergedBuf = NULL;
for (block = blockList; block != NULL; )
    {
    /* Find contigious blocks and read them into mergedBuf. */
    fileOffsetSizeFindGap(block, &beforeGap, &afterGap);
    bits64 mergedOffset = block->offset;
    bits64 mergedSize = beforeGap->offset + beforeGap->size - mergedOffset;
    udcSeek(udc, mergedOffset);
    mergedBuf = needLargeMem(mergedSize);
    udcMustRead(udc, mergedBuf, mergedSize);
    char *blockBuf = mergedBuf;

    /* Loop through individual blocks within merged section. */
    for (;block != afterGap; block = block->next)
        {
	/* Uncompress if necessary. */
	char *blockPt, *blockEnd;
	if (uncompressBuf)
	    {
	    blockPt = uncompressBuf;
	    int uncSize = zUncompress(blockBuf, block->size, uncompressBuf, bbi->uncompressBufSize);
	    blockEnd = blockPt + uncSize;
	    }
	else
	    {
	    blockPt = blockBuf;
	    blockEnd = blockPt + block->size;
	    }

	while (blockPt < blockEnd)
	    {
	    /* Read next record into local variables. */
	    bits32 chr = memReadBits32(&blockPt, isSwapped);
	    bits32 s = memReadBits32(&blockPt, isSwapped);
	    bits32 e = memReadBits32(&blockPt, isSwapped);

	    /* calculate length of rest of bed fields */
	    int restLen = strlen(blockPt);

	    /* If we're actually in range then copy it into a new  element and add to list. */
	    if (chr == chromId &&
                ((s < end && e > start)
                // Make sure to include zero-length insertion elements at start or end:
                 || (s == e && (s == end || e == start))))
		{
		++itemCount;
		if (maxItems > 0 && itemCount > maxItems)
		    break;

		lmAllocVar(lm, el);
		el->start = s;
		el->end = e;
		if (restLen > 0)
		    el->rest = lmCloneStringZ(lm, blockPt, restLen);
		el->chromId = chromId;
		slAddHead(&list, el);
		}

	    // move blockPt pointer to end of previous bed
	    blockPt += restLen + 1;
	    }
	if (maxItems > 0 && itemCount > maxItems)
	    break;
	blockBuf += block->size;
        }
    if (maxItems > 0 && itemCount > maxItems)
        break;
    freez(&mergedBuf);
    }
freez(&mergedBuf);
freeMem(uncompressBuf);
slFreeList(&blockList);
slReverse(&list);
return list;
}

int bigBedIntervalToRow(struct bigBedInterval *interval, char *chrom, char *startBuf, char *endBuf,
	char **row, int rowSize)
/* Convert bigBedInterval into an array of chars equivalent to what you'd get by
 * parsing the bed file. The startBuf and endBuf are used to hold the ascii representation of
 * start and end.  Note that the interval->rest string will have zeroes inserted as a side effect.
 */
{
int fieldCount = 3;
snprintf(startBuf, 10, "%u", interval->start);
snprintf(endBuf, 10, "%u", interval->end);
row[0] = chrom;
row[1] = startBuf;
row[2] = endBuf;
if (!isEmpty(interval->rest))
    {
    int wordCount = chopByChar(cloneString(interval->rest), '\t', row+3, rowSize-3);
    fieldCount += wordCount;
    }
return fieldCount;
}

struct offsetSize 
/* Simple file offset and file size. */
    {
    bits64 offset; 
    bits64 size;
    };

typedef boolean (*BbFirstWordMatch)(char *line, int fieldIx, void *target);
/* A function that returns TRUE if first word in tab-separated line matches target. */

char *bigBedAutoSqlText(struct bbiFile *bbi)
/* Get autoSql text if any associated with file.  Do a freeMem of this when done. */
{
if (bbi->asOffset == 0)
    return NULL;
struct udcFile *f = bbi->udc;
udcSeek(f, bbi->asOffset);
return udcReadStringAndZero(f);
}

struct asObject *bigBedAs(struct bbiFile *bbi)
/* Get autoSql object definition if any associated with file. */
{
if (bbi->asOffset == 0)
    return NULL;
char *asText = bigBedAutoSqlText(bbi);
struct asObject *as = asParseText(asText);
freeMem(asText);
return as;
}

struct asObject *bigBedAsOrDefault(struct bbiFile *bbi)
// Get asObject associated with bigBed - if none exists in file make it up from field counts.
{
struct asObject *as = bigBedAs(bbi);
if (as == NULL)
    as = asParseText(bedAsDef(bbi->definedFieldCount, bbi->fieldCount));
return as;
}
