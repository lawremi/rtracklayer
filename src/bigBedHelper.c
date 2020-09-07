#include "bigBedHelper.h"

/*
  Most of the functions in this file are taken from ucscGenomeBrowser/kent/src/utils/bedToBigBed.c
*/

int getDefinedFieldCount(struct asObject *as) {
  int definedFieldCount = 0;
  struct asColumn *asCol = as->columnList;
  char *asText = bedAsDef(12, 12);
  struct asObject *bedAs = asParseText(asText);
  freeMem(asText);
  struct asColumn *bedCol = bedAs->columnList;
  while (asCol && bedCol) {
    if (strncmp(asCol->name, bedCol->name, strlen(asCol->name)) == 0)
      ++definedFieldCount;
    bedCol = bedCol->next;
    asCol = asCol->next;
  }
  asObjectFree(&bedAs);
  return definedFieldCount;
}

bool isPresent(int definedFieldCount, int index) {
  return (definedFieldCount - index >= 0) ? TRUE : FALSE;
}

bool isSelected(SEXP r_selectedindex, int position) {
  if (length(r_selectedindex) == 0)
    return TRUE;
  for (int i = 0; i < length(r_selectedindex); ++i) {
    if(INTEGER(r_selectedindex)[i] == position)
      return TRUE;
  }
  return FALSE;
}

/* following functions are taken from the kent library */

struct rbTree *rangeTreeForBedChrom(struct lineFile *lf, char *chrom)
/* Read lines from bed file as long as they match chrom.  Return a rangeTree that
 * corresponds to the coverage. */
{
struct rbTree *tree = rangeTreeNew();
char *line;
while (lineFileNextReal(lf, &line))
    {
    if (!startsWithWord(chrom, line))
        {
	lineFileReuse(lf);
	break;
	}
    char *row[3];
    chopLine(line, row);
    unsigned start = sqlUnsigned(row[1]);
    unsigned end = sqlUnsigned(row[2]);
    rangeTreeAddToCoverageDepth(tree, start, end);
    }
return tree;
}

void bbExIndexMakerAddKeysFromRow(struct bbExIndexMaker *eim, char **row, int recordIx)
/* Save the keys that are being indexed by row in eim. */
{
int i;
for (i=0; i < eim->indexCount; ++i)
    {
    int rowIx = eim->indexFields[i];
    eim->chunkArrayArray[i][recordIx].name = cloneString(row[rowIx]);
    }
}

void bbExIndexMakerAddOffsetSize(struct bbExIndexMaker *eim, bits64 offset, bits64 size,
    long startIx, long endIx)
/* Update offset and size fields of all file chunks between startIx and endIx */
{
int i;
for (i=0; i < eim->indexCount; ++i)
    {
    struct bbNamedFileChunk *chunks = eim->chunkArrayArray[i];
    long j;
    for (j = startIx; j < endIx; ++j)
        {
	struct bbNamedFileChunk *chunk = chunks + j;
	chunk->offset = offset;
	chunk->size = size;
	}
    }
}

/* Allocate the big part of the extra index maker - the part that holds which
 * chunk is used for each record. */
void bbExIndexMakerAllocChunkArrays(struct bbExIndexMaker *eim, int recordCount) {
  eim->recordCount = recordCount;
  int i;
  for (i=0; i < eim->indexCount; ++i)
    AllocArray(eim->chunkArrayArray[i], recordCount);
}

/* Return an index maker corresponding to extraIndexList. Checks that all fields
 * mentioned are in autoSql definition, and for now that they are all text fields. */
struct bbExIndexMaker *bbExIndexMakerNew(struct slName *extraIndexList, struct asObject *as) {
  /* Fill in scalar fields and return quickly if no extra indexes. */
  struct bbExIndexMaker *eim;
  AllocVar(eim);
  eim->indexCount = slCount(extraIndexList);
  if (eim->indexCount == 0)
    return eim;	// Not much to do in this case

  /* Allocate arrays according field count. */
  AllocArray(eim->indexFields, eim->indexCount);
  AllocArray(eim->maxFieldSize, eim->indexCount);
  AllocArray(eim->chunkArrayArray, eim->indexCount);
  AllocArray(eim->fileOffsets, eim->indexCount);

  /* Loop through each field checking that it is indeed something we can index
  * and if so saving information about it */
  int indexIx = 0;
  struct slName *name;
  for (name = extraIndexList; name != NULL; name = name->next) {
      struct asColumn *col = asColumnFind(as, name->name);
      if (col == NULL)
          errAbort("extraIndex field %s not a standard bed field or found in autoSql string.",
        name->name);
      if (!sameString(col->lowType->name, "string"))
          errAbort("Sorry for now can only index string fields.");
      eim->indexFields[indexIx] = slIxFromElement(as->columnList, col);
      ++indexIx;
  }
  return eim;
}

/* Compare two named offset object to facilitate qsorting by name. */
int bbNamedFileChunkCmpByName(const void *va, const void *vb) {
  const struct bbNamedFileChunk *a = va, *b = vb;
  return strcmp(a->name, b->name);
}

/* Return pointer to val for bPlusTree maker. */
void *bbNamedFileChunkVal(const void *va) {
  const struct bbNamedFileChunk *item = va;
  return (void *)&item->offset;
}

/* Copy name to keyBuf for bPlusTree maker */
void bbNamedFileChunkKey(const void *va, char *keyBuf) {
  const struct bbNamedFileChunk *item = va;
  strncpy(keyBuf,item->name, maxBedNameSize);
}

struct bbiSummary *bedWriteReducedOnceReturnReducedTwice(struct bbiChromUsage *usageList,
	int fieldCount, struct lineFile *lf, bits32 initialReduction, bits32 initialReductionCount,
	int zoomIncrement, int blockSize, int itemsPerSlot, boolean doCompress,
	struct lm *lm, FILE *f, bits64 *retDataStart, bits64 *retIndexStart,
	struct bbiSummaryElement *totalSum)
/* Write out data reduced by factor of initialReduction.  Also calculate and keep in memory
 * next reduction level.  This is more work than some ways, but it keeps us from having to
 * keep the first reduction entirely in memory. */
{
struct bbiSummary *twiceReducedList = NULL;
bits32 doubleReductionSize = initialReduction * zoomIncrement;
struct bbiChromUsage *usage = usageList;
struct bbiBoundsArray *boundsArray, *boundsPt, *boundsEnd;
boundsPt = AllocArray(boundsArray, initialReductionCount);
boundsEnd = boundsPt + initialReductionCount;

*retDataStart = ftell(f);
writeOne(f, initialReductionCount);
/* This gets a little complicated I'm afraid.  The strategy is to:
 *   1) Build up a range tree that represents coverage depth on that chromosome
 *      This also has the nice side effect of getting rid of overlaps.
 *   2) Stream through the range tree, outputting the initial summary level and
 *      further reducing. 
 */
boolean firstTime = TRUE;
struct bbiSumOutStream *stream = bbiSumOutStreamOpen(itemsPerSlot, f, doCompress);
for (usage = usageList; usage != NULL; usage = usage->next)
    {
    struct bbiSummary oneSummary, *sum = NULL;
    struct rbTree *rangeTree = rangeTreeForBedChrom(lf, usage->name);
    struct range *range, *rangeList = rangeTreeList(rangeTree);
    for (range = rangeList; range != NULL; range = range->next)
        {
	/* Grab values we want from range. */
	double val = ptToInt(range->val);
	int start = range->start;
	int end = range->end;
	bits32 size = end - start;

	/* Add to total summary. */
	if (firstTime)
	    {
	    totalSum->validCount = size;
	    totalSum->minVal = totalSum->maxVal = val;
	    totalSum->sumData = val*size;
	    totalSum->sumSquares = val*val*size;
	    firstTime = FALSE;
	    }
	else
	    {
	    totalSum->validCount += size;
	    if (val < totalSum->minVal) totalSum->minVal = val;
	    if (val > totalSum->maxVal) totalSum->maxVal = val;
	    totalSum->sumData += val*size;
	    totalSum->sumSquares += val*val*size;
	    }

	/* If start past existing block then output it. */
	if (sum != NULL && sum->end <= start && sum->end < usage->size)
	    {
	    bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize,
		&boundsPt, boundsEnd, lm, stream);
	    sum = NULL;
	    }
	/* If don't have a summary we're working on now, make one. */
	if (sum == NULL)
	    {
	    oneSummary.chromId = usage->id;
	    oneSummary.start = start;
	    oneSummary.end = start + initialReduction;
	    if (oneSummary.end > usage->size) oneSummary.end = usage->size;
	    oneSummary.minVal = oneSummary.maxVal = val;
	    oneSummary.sumData = oneSummary.sumSquares = 0.0;
	    oneSummary.validCount = 0;
	    sum = &oneSummary;
	    }
	/* Deal with case where might have to split an item between multiple summaries.  This
	 * loop handles all but the final affected summary in that case. */
	while (end > sum->end)
	    {
	    /* Fold in bits that overlap with existing summary and output. */
	    int overlap = rangeIntersection(start, end, sum->start, sum->end);
	    assert(overlap > 0);
	    sum->validCount += overlap;
	    if (sum->minVal > val) sum->minVal = val;
	    if (sum->maxVal < val) sum->maxVal = val;
	    sum->sumData += val * overlap;
	    sum->sumSquares += val*val * overlap;
	    bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize,
		    &boundsPt, boundsEnd, lm, stream);
	    size -= overlap;

	    /* Move summary to next part. */
	    sum->start = start = sum->end;
	    sum->end = start + initialReduction;
	    if (sum->end > usage->size) sum->end = usage->size;
	    sum->minVal = sum->maxVal = val;
	    sum->sumData = sum->sumSquares = 0.0;
	    sum->validCount = 0;
	    }

	/* Add to summary. */
	sum->validCount += size;
	if (sum->minVal > val) sum->minVal = val;
	if (sum->maxVal < val) sum->maxVal = val;
	sum->sumData += val * size;
	sum->sumSquares += val*val * size;
	}
    if (sum != NULL)
	{
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize,
	    &boundsPt, boundsEnd, lm, stream);
	}
    rangeTreeFree(&rangeTree);
    }
bbiSumOutStreamClose(&stream);

/* Write out 1st zoom index. */
int indexOffset = *retIndexStart = ftell(f);
assert(boundsPt == boundsEnd);
cirTreeFileBulkIndexToOpenFile(boundsArray, sizeof(boundsArray[0]), initialReductionCount,
    blockSize, itemsPerSlot, NULL, bbiBoundsArrayFetchKey, bbiBoundsArrayFetchOffset,
    indexOffset, f);

freez(&boundsArray);
slReverse(&twiceReducedList);
return twiceReducedList;
}

void writeBlocks(struct bbiChromUsage *usageList, struct lineFile *lf, struct asObject *as,
	int itemsPerSlot, struct bbiBoundsArray *bounds,
	int sectionCount, boolean doCompress, FILE *f,
	int resTryCount, int resScales[], int resSizes[],
	struct bbExIndexMaker *eim,  int bedCount,
	bits16 fieldCount, int bedN, bits32 *retMaxBlockSize)
/* Read through lf, writing it in f.  Save starting points of blocks (every itemsPerSlot)
 * to boundsArray */
{
int maxBlockSize = 0;
struct bbiChromUsage *usage = usageList;
char *line, *row[fieldCount+1];
int lastField = fieldCount-1;
int itemIx = 0, sectionIx = 0;
bits64 blockStartOffset = 0;
int startPos = 0, endPos = 0;
bits32 chromId = 0;
struct dyString *stream = dyStringNew(0);

/* Will keep track of some things that help us determine how much to reduce. */
bits32 resEnds[resTryCount];
int resTry;
for (resTry = 0; resTry < resTryCount; ++resTry)
    resEnds[resTry] = 0;
boolean atEnd = FALSE, sameChrom = FALSE;
bits32 start = 0, end = 0;
char *chrom = NULL;
struct bed *bed;
AllocVar(bed);

/* Help keep track of which beds are in current chunk so as to write out
 * namedChunks to eim if need be. */
long sectionStartIx = 0, sectionEndIx = 0;

for (;;)
    {
    /* Get next line of input if any. */
    if (lineFileNextReal(lf, &line))
	{
	/* Chop up line and make sure the word count is right. */
	int wordCount;
	wordCount = chopLine(line, row);
	lineFileExpectWords(lf, fieldCount, wordCount);
	loadAndValidateBedExt(row, bedN, fieldCount, lf, bed, as, FALSE, TRUE);
	chrom = bed->chrom;
	start = bed->chromStart;
	end = bed->chromEnd;

	sameChrom = sameString(chrom, usage->name);
	}
    else  /* No next line */
	{
	atEnd = TRUE;
	}

    /* Check conditions that would end block and save block info and advance to next if need be. */
    if (atEnd || !sameChrom || itemIx >= itemsPerSlot)
        {
	/* Save stream to file, compressing if need be. */
	if (stream->stringSize > maxBlockSize)
	    maxBlockSize = stream->stringSize;
	if (doCompress)
            {
	    size_t maxCompSize = zCompBufSize(stream->stringSize);

            // keep around an area of scratch memory
            static int compBufSize = 0;
            static char *compBuf = NULL;
            // check to see if buffer needed for compression is big enough
            if (compBufSize < maxCompSize)
                {
                // free up the old not-big-enough piece
                freez(&compBuf); // freez knows bout NULL

                // get new scratch area
                compBufSize = maxCompSize;
                compBuf = needLargeMem(compBufSize);
                }
	    int compSize = zCompress(stream->string, stream->stringSize, compBuf, maxCompSize);
	    mustWrite(f, compBuf, compSize);
	    }
	else
	    mustWrite(f, stream->string, stream->stringSize);
	dyStringClear(stream);

	/* Save block offset and size for all named chunks in this section. */
	if (eim != NULL)
	    {
	    bits64 blockEndOffset = ftell(f);
	    bbExIndexMakerAddOffsetSize(eim, blockStartOffset, blockEndOffset-blockStartOffset,
		sectionStartIx, sectionEndIx);
	    sectionStartIx = sectionEndIx;
	    }

	/* Save info on existing block. */
	struct bbiBoundsArray *b = &bounds[sectionIx];
	b->offset = blockStartOffset;
	b->range.chromIx = chromId;
	b->range.start = startPos;
	b->range.end = endPos;
	++sectionIx;
	itemIx = 0;

	if (atEnd)
	    break;
	}

    /* Advance to next chromosome if need be and get chromosome id. */
    if (!sameChrom)
        {
	usage = usage->next;
	assert(usage != NULL);
	assert(sameString(chrom, usage->name));
	for (resTry = 0; resTry < resTryCount; ++resTry)
	    resEnds[resTry] = 0;
	}
    chromId = usage->id;

    /* At start of block we save a lot of info. */
    if (itemIx == 0)
        {
	blockStartOffset = ftell(f);
	startPos = start;
	endPos = end;
	}
    /* Otherwise just update end. */
        {
	if (endPos < end)
	    endPos = end;
	/* No need to update startPos since list is sorted. */
	}

    /* Save name into namedOffset if need be. */
    if (eim != NULL)
	{
	bbExIndexMakerAddKeysFromRow(eim, row, sectionEndIx);
	sectionEndIx += 1;
	}

    /* Write out data. */
    dyStringWriteOne(stream, chromId);
    dyStringWriteOne(stream, start);
    dyStringWriteOne(stream, end);
    if (fieldCount > 3)
        {
	int i;
	/* Write 3rd through next to last field and a tab separator. */
	for (i=3; i<lastField; ++i)
	    {
	    char *s = row[i];
	    dyStringAppend(stream, s);
	    dyStringAppendC(stream, '\t');
	    }
	/* Write last field and terminal zero */
	char *s = row[lastField];
	dyStringAppend(stream, s);
	}
    dyStringAppendC(stream, 0);

    itemIx += 1;

    /* Do zoom counting. */
    for (resTry = 0; resTry < resTryCount; ++resTry)
        {
	bits32 resEnd = resEnds[resTry];
	if (start >= resEnd)
	    {
	    resSizes[resTry] += 1;
	    resEnds[resTry] = resEnd = start + resScales[resTry];
	    }
	while (end > resEnd)
	    {
	    resSizes[resTry] += 1;
	    resEnds[resTry] = resEnd = resEnd + resScales[resTry];
	    }
	}
    }
assert(sectionIx == sectionCount);
freez(&bed);
*retMaxBlockSize = maxBlockSize;
}
