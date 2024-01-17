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
sprintf(startBuf, "%u", interval->start);
sprintf(endBuf, "%u", interval->end);
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

static struct bbiInterval *bigBedCoverageIntervals(struct bbiFile *bbi,
	char *chrom, bits32 start, bits32 end, struct lm *lm)
/* Return intervals where the val is the depth of coverage. */
{
/* Get list of overlapping intervals */
struct bigBedInterval *bi, *biList = bigBedIntervalQuery(bbi, chrom, start, end, 0, lm);
if (biList == NULL)
    return NULL;

/* Make a range tree that collects coverage. */
struct rbTree *rangeTree = rangeTreeNew();
for (bi = biList; bi != NULL; bi = bi->next)
    rangeTreeAddToCoverageDepth(rangeTree, bi->start, bi->end);
struct range *range, *rangeList = rangeTreeList(rangeTree);

/* Convert rangeList to bbiInterval list. */
struct bbiInterval *bwi, *bwiList = NULL;
for (range = rangeList; range != NULL; range = range->next)
    {
    lmAllocVar(lm, bwi);
    bwi->start = range->start;
    if (bwi->start < start)
       bwi->start = start;
    bwi->end = range->end;
    if (bwi->end > end)
       bwi->end = end;
    bwi->val = ptToInt(range->val);
    slAddHead(&bwiList, bwi);
    }
slReverse(&bwiList);

/* Clean up and go home. */
rangeTreeFree(&rangeTree);
return bwiList;
}


struct offsetSize 
/* Simple file offset and file size. */
    {
    bits64 offset; 
    bits64 size;
    };

static int cmpOffsetSizeRef(const void *va, const void *vb)
/* Compare to sort slRef pointing to offsetSize.  Sort is kind of hokey,
 * but guarantees all items that are the same will be next to each other
 * at least, which is all we care about. */
{
const struct slRef *a = *((struct slRef **)va);
const struct slRef *b = *((struct slRef **)vb);
return memcmp(a->val, b->val, sizeof(struct offsetSize));
}

static struct fileOffsetSize *fosFromRedundantBlockList(struct slRef **pBlockList, 
    boolean isSwapped)
/* Convert from list of references to offsetSize format to list of fileOffsetSize
 * format, while removing redundancy.   Sorts *pBlockList as a side effect. */
{
/* Sort input so it it easy to uniquify. */
slSort(pBlockList, cmpOffsetSizeRef);
struct slRef *blockList = *pBlockList;

/* Make new fileOffsetSize for each unique offsetSize. */
struct fileOffsetSize *fosList = NULL, *fos;
struct offsetSize lastOffsetSize = {0,0};
struct slRef *blockRef;
for (blockRef = blockList; blockRef != NULL; blockRef = blockRef->next)
    {
    if (memcmp(&lastOffsetSize, blockRef->val, sizeof(lastOffsetSize)) != 0)
        {
	memcpy(&lastOffsetSize, blockRef->val, sizeof(lastOffsetSize));
	AllocVar(fos);
	if (isSwapped)
	    {
	    fos->offset = byteSwap64(lastOffsetSize.offset);
	    fos->size = byteSwap64(lastOffsetSize.size);
	    }
	else
	    {
	    fos->offset = lastOffsetSize.offset;
	    fos->size = lastOffsetSize.size;
	    }
	slAddHead(&fosList, fos);
	}
    }
slReverse(&fosList);
return fosList;
}


static struct fileOffsetSize *bigBedChunksMatchingName(struct bbiFile *bbi, 
    struct bptFile *index, char *name)
/* Get list of file chunks that match name.  Can slFreeList this when done. */
{
struct slRef *blockList = bptFileFindMultiple(index, 
	name, strlen(name), sizeof(struct offsetSize));
struct fileOffsetSize *fosList = fosFromRedundantBlockList(&blockList, bbi->isSwapped);
slRefFreeListAndVals(&blockList);
return fosList;
}

static struct fileOffsetSize *bigBedChunksMatchingNames(struct bbiFile *bbi, 
	struct bptFile *index, char **names, int nameCount)
/* Get list of file chunks that match any of the names.  Can slFreeList this when done. */
{
/* Go through all names and make a blockList that includes all blocks with any hit to any name.  
 * Many of these blocks will occur multiple times. */
struct slRef *blockList = NULL;
int nameIx;
for (nameIx = 0; nameIx < nameCount; ++nameIx)
    {
    char *name = names[nameIx];
    struct slRef *oneList = bptFileFindMultiple(index, 
	    name, strlen(name), sizeof(struct offsetSize));
    blockList = slCat(oneList, blockList);
    }

/* Create nonredundant list of blocks. */
struct fileOffsetSize *fosList = fosFromRedundantBlockList(&blockList, bbi->isSwapped);

/* Clean up and resturn result. */
slRefFreeListAndVals(&blockList);
return fosList;
}

typedef boolean (*BbFirstWordMatch)(char *line, int fieldIx, void *target);
/* A function that returns TRUE if first word in tab-separated line matches target. */

static void extractField(char *line, int fieldIx, char **retField, int *retFieldSize)
/* Go through tab separated line and figure out start and size of given field. */
{
int i;
fieldIx -= 3;	/* Skip over chrom/start/end, which are not in line. */
for (i=0; i<fieldIx; ++i)
    {
    line = strchr(line, '\t');
    if (line == NULL)
        {
	warn("Not enough fields in extractField of %s", line);
	internalErr();
	}
    line += 1;
    }
char *end = strchr(line, '\t');
if (end == NULL)
    end = line + strlen(line);
*retField = line;
*retFieldSize = end - line;
}

static boolean bbWordMatchesName(char *line, int fieldIx, void *target)
/* Return true if first word of line is same as target, which is just a string. */
{
char *name = target;
int fieldSize;
char *field;
extractField(line, fieldIx, &field, &fieldSize);
return strlen(name) == fieldSize && memcmp(name, field, fieldSize) == 0;
}

static boolean bbWordIsInHash(char *line, int fieldIx, void *target)
/* Return true if first word of line is same as target, which is just a string. */
{
int fieldSize;
char *field;
extractField(line, fieldIx, &field, &fieldSize);
char fieldString[fieldSize+1];
memcpy(fieldString, field, fieldSize);
fieldString[fieldSize] = 0;

/* Return boolean value that reflects whether we found it in hash */
struct hash *hash = target;
return hashLookup(hash, fieldString) != NULL;
}

static struct bigBedInterval *bigBedIntervalsMatchingName(struct bbiFile *bbi, 
    struct fileOffsetSize *fosList, BbFirstWordMatch matcher, int fieldIx, 
    void *target, struct lm *lm)
/* Return list of intervals inside of sectors of bbiFile defined by fosList where the name 
 * matches target somehow. */
{
struct bigBedInterval *interval, *intervalList = NULL;
struct fileOffsetSize *fos;
boolean isSwapped = bbi->isSwapped;
for (fos = fosList; fos != NULL; fos = fos->next)
    {
    /* Read in raw data */
    udcSeek(bbi->udc, fos->offset);
    char *rawData = needLargeMem(fos->size);
    udcRead(bbi->udc, rawData, fos->size);

    /* Optionally uncompress data, and set data pointer to uncompressed version. */
    char *uncompressedData = NULL;
    char *data = NULL;
    int dataSize = 0;
    if (bbi->uncompressBufSize > 0)
	{
	data = uncompressedData = needLargeMem(bbi->uncompressBufSize);
	dataSize = zUncompress(rawData, fos->size, uncompressedData, bbi->uncompressBufSize);
	}
    else
	{
        data = rawData;
	dataSize = fos->size;
	}

    /* Set up for "memRead" routines to more or less treat memory block like file */
    char *blockPt = data, *blockEnd = data + dataSize;
    struct dyString *dy = dyStringNew(32); // Keep bits outside of chrom/start/end here


    /* Read next record into local variables. */
    while (blockPt < blockEnd)
	{
	bits32 chromIx = memReadBits32(&blockPt, isSwapped);
	bits32 s = memReadBits32(&blockPt, isSwapped);
	bits32 e = memReadBits32(&blockPt, isSwapped);
	int c;
	dyStringClear(dy);
	// TODO - can simplify this probably just to for (;;) {if ((c = *blockPt++) == 0) ...
	while ((c = *blockPt++) >= 0)
	    {
	    if (c == 0)
		break;
	    dyStringAppendC(dy, c);
	    }
	if ((*matcher)(dy->string, fieldIx, target))
	    {
	    lmAllocVar(lm, interval);
	    interval->start = s;
	    interval->end = e;
	    interval->rest = cloneString(dy->string);
	    interval->chromId = chromIx;
	    slAddHead(&intervalList, interval);
	    }
	}

    /* Clean up temporary buffers. */
    dyStringFree(&dy);
    freez(&uncompressedData);
    freez(&rawData);
    }
slReverse(&intervalList);
return intervalList;
}


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
