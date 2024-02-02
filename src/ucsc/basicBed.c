/* basicBed contains the basic code for Browser Extensible Data (bed) files and tables.
 * The idea behind bed is that the first three fields are defined and required.
 * A total of 15 fields are defined, and the file can contain any number of these.
 * In addition after any number of defined fields there can be custom fields that
 * are not defined in the bed spec.
 *
 * There's additional bed-related code in src/hg/inc/bed.h.  This module contains the
 * stuff that's independent of the database and other genomic structures. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */


#include "common.h"
#include "hash.h"
#include "linefile.h"
#include "dystring.h"
#include "sqlNum.h"
#include "sqlList.h"
#include "rangeTree.h"
#include "binRange.h"
#include "asParse.h"
#include "htmlColor.h"
#include "basicBed.h"
#include "memgfx.h"
#include "localmem.h"

void bedFree(struct bed **pEl)
/* Free a single dynamically allocated bed such as created
 * with bedLoad(). */
{
struct bed *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->blockSizes);
freeMem(el->chromStarts);
freeMem(el->expIds);
freeMem(el->expScores);
freez(pEl);
}

/* --------------- End of AutoSQL generated code. --------------- */


struct bedLine *bedLineNew(char *line)
/* Create a new bedLine based on tab-separated string s. */
{
struct bedLine *bl;
char *s, c;

AllocVar(bl);
bl->chrom = cloneString(line);
s = strchr(bl->chrom, '\t');
if (s == NULL)
    errAbort("Expecting tab in bed line %s", line);
*s++ = 0;
c = *s;
if (isdigit(c) || (c == '-' && isdigit(s[1])))
    {
    bl->chromStart = atoi(s);
    bl->line = s;
    return bl;
    }
else
    {
    errAbort("Expecting start position in second field of %s", line);
    return NULL;
    }
}

void bedLineFree(struct bedLine **pBl)
/* Free up memory associated with bedLine. */
{
struct bedLine *bl;

if ((bl = *pBl) != NULL)
    {
    freeMem(bl->chrom);
    freez(pBl);
    }
}

int bedLineCmp(const void *va, const void *vb)
/* Compare to sort based on query. */
{
const struct bedLine *a = *((struct bedLine **)va);
const struct bedLine *b = *((struct bedLine **)vb);
int dif;
dif = strcmp(a->chrom, b->chrom);
if (dif == 0)
    dif = a->chromStart - b->chromStart;
return dif;
}

struct bed *bedLoad5(char **row)
/* Load first five fields of bed. */
{
struct bed *ret;
AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlSigned(row[4]);
return ret;
}

/* it turns out that it isn't just hgLoadBed and custom tracks
 *	that may encounter the r,g,b specification.  Any program that
 *	reads bed files may enconter them, so take care of them
 *	at any time.  The strchr() function is very fast which will
 *	be a failure in the vast majority of cases parsing integers,
 *	therefore, this shouldn't be too severe a performance hit.
 */
int itemRgbColumn(char *column9)
/* Convert color specification to internal format. */
{
int itemRgb = 0;
/*  Allow comma separated list of rgb values here   */
char *comma = strchr(column9, ',');
if (comma)
    {
    if (-1 == (itemRgb = bedParseRgb(column9)))
	errAbort("ERROR: expecting r,g,b specification, "
		    "found: '%s'", column9);
    }
else
    itemRgb = sqlUnsigned(column9);
return itemRgb;
}

struct bed *bedLoadN(char *row[], int wordCount)
/* Convert a row of strings to a bed. */
{
struct bed * bed;
int count;

AllocVar(bed);
bed->chrom = cloneString(row[0]);
bed->chromStart = sqlUnsigned(row[1]);
bed->chromEnd = sqlUnsigned(row[2]);
if (wordCount > 3)
     bed->name = cloneString(row[3]);
if (wordCount > 4)
     bed->score = sqlSigned(row[4]);
if (wordCount > 5)
     bed->strand[0] = row[5][0];
if (wordCount > 6)
     bed->thickStart = sqlUnsigned(row[6]);
else
     bed->thickStart = bed->chromStart;
if (wordCount > 7)
     bed->thickEnd = sqlUnsigned(row[7]);
else
     bed->thickEnd = bed->chromEnd;
if (wordCount > 8)
    bed->itemRgb = itemRgbColumn(row[8]);
if (wordCount > 9)
    bed->blockCount = sqlUnsigned(row[9]);
if (wordCount > 10)
    sqlSignedDynamicArray(row[10], &bed->blockSizes, &count);
if (wordCount > 11)
    sqlSignedDynamicArray(row[11], &bed->chromStarts, &count);
if (wordCount > 12)
    bed->expCount = sqlUnsigned(row[12]);
if (wordCount > 13)
    sqlSignedDynamicArray(row[13], &bed->expIds, &count);
if (wordCount > 14)
    sqlFloatDynamicArray(row[14], &bed->expScores, &count);
return bed;
}

struct bed *bedLoadNAllChrom(char *fileName, int numFields, char* chrom) 
/* Load bed entries from a tab-separated file that have the given chrom.
 * Dispose of this with bedFreeList(). */
{
struct bed *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[numFields];

while (lineFileRow(lf, row))
    {
    el = bedLoadN(row, numFields);
    if(chrom == NULL || sameString(el->chrom, chrom))
        slAddHead(&list, el);
    else
        bedFree(&el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}


void bedLoadAllReturnFieldCountAndRgb(char *fileName, struct bed **retList, int *retFieldCount, 
    boolean *retRgb)
/* Load bed of unknown size and return number of fields as well as list of bed items.
 * Ensures that all lines in bed file have same field count.  Also returns whether 
 * column 9 is being used as RGB or not. */
{
struct bed *list = NULL;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *line, *row[bedKnownFields];
int fieldCount = 0;
boolean isRgb = FALSE;

while (lineFileNextReal(lf, &line))
    {
    int numFields = chopByWhite(line, row, ArraySize(row));
    if (numFields < 4)
	errAbort("file %s doesn't appear to be in bed format. At least 4 fields required, got %d", 
		fileName, numFields);
    if (fieldCount == 0)
	{
        fieldCount = numFields;
	if (fieldCount >= 9)
	    isRgb =  (strchr(row[8], ',') != NULL);
	}
    else
        if (fieldCount != numFields)
	    errAbort("Inconsistent number of fields in file. %d on line %d of %s, %d previously.",
	        numFields, lf->lineIx, lf->fileName, fieldCount);
    slAddHead(&list, bedLoadN(row, fieldCount));
    }
lineFileClose(&lf);
slReverse(&list);
*retList = list;
*retFieldCount = fieldCount;
if (retRgb != NULL)
   *retRgb = isRgb;
}


void bedOutFlexible(struct bed *el, int wordCount, FILE *f,
	char sep, char lastSep, boolean useItemRgb)
/* Write a bed of wordCount fields, optionally interpreting field nine as R,G,B values. */
{
int i;
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->chromStart);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
if (wordCount <= 3)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
if (wordCount <= 4)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
fprintf(f, "%d", el->score);
if (wordCount <= 5)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->strand);
if (sep == ',') fputc('"',f);
if (wordCount <= 6)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
fprintf(f, "%u", el->thickStart);
if (wordCount <= 7)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
fprintf(f, "%u", el->thickEnd);
if (wordCount <= 8)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (useItemRgb)
    fprintf(f, "%d,%d,%d", (el->itemRgb & 0xff0000) >> 16,
        (el->itemRgb & 0xff00) >> 8, (el->itemRgb & 0xff));
else
    fprintf(f, "%u", el->itemRgb);
if (wordCount <= 9)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
fprintf(f, "%d", el->blockCount);
if (wordCount <= 10)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->blockCount; ++i)
    {
    fprintf(f, "%d", el->blockSizes[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
if (wordCount <= 11)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->blockCount; ++i)
    {
    fprintf(f, "%d", el->chromStarts[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);

if (wordCount <= 12)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
fprintf(f, "%d", el->expCount);

if (wordCount <= 13)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->expCount; ++i)
    {
    fprintf(f, "%d", el->expIds[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);


if (wordCount <= 14)
    {
    fputc(lastSep, f);
    return;
    }
fputc(sep,f);
if (sep == ',') fputc('{',f);
for (i=0; i<el->expCount; ++i)
    {
    fprintf(f, "%g", el->expScores[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);


fputc(lastSep,f);
}


int bedTotalBlockSize(struct bed *bed)
/* Return total size of all blocks. */
{
int total = 0;
int i;
if (bed->blockCount == 0)
    return bed->chromEnd - bed->chromStart;
for (i=0; i<bed->blockCount; ++i)
    total += bed->blockSizes[i];
return total;
}

int bedBlockSizeInRange(struct bed *bed, int rangeStart, int rangeEnd)
/* Get size of all parts of all exons between rangeStart and rangeEnd. */
{
int total = 0;
int i;
for (i=0; i<bed->blockCount; ++i)
    {
    int start = bed->chromStart + bed->chromStarts[i];
    int end = start + bed->blockSizes[i];
    total += positiveRangeIntersection(start, end, rangeStart, rangeEnd);
    }
return total;
}


struct bed *cloneBed(struct bed *bed)
/* Make an all-newly-allocated copy of a single bed record. */
{
struct bed *newBed;
if (bed == NULL)
    return NULL;
AllocVar(newBed);
newBed->chrom = cloneString(bed->chrom);
newBed->chromStart = bed->chromStart;
newBed->chromEnd = bed->chromEnd;
newBed->name = cloneString(bed->name);
newBed->score = bed->score;
strncpy(newBed->strand, bed->strand, sizeof(bed->strand));
newBed->thickStart = bed->thickStart;
newBed->thickEnd = bed->thickEnd;
newBed->itemRgb = bed->itemRgb;
newBed->blockCount = bed->blockCount;
if (bed->blockCount > 0)
    {
    newBed->blockSizes = needMem(sizeof(int) * bed->blockCount);
    memcpy(newBed->blockSizes, bed->blockSizes,
	   sizeof(int) * bed->blockCount);
    newBed->chromStarts = needMem(sizeof(int) * bed->blockCount);
    memcpy(newBed->chromStarts, bed->chromStarts,
	   sizeof(int) * bed->blockCount);
    }
newBed->expCount = bed->expCount;
if (bed->expCount > 0)
    {
    newBed->expIds = needMem(sizeof(int) * bed->expCount);
    memcpy(newBed->expIds, bed->expIds,
	   sizeof(int) * bed->expCount);
    newBed->expScores = needMem(sizeof(float) * bed->expCount);
    memcpy(newBed->expScores, bed->expScores,
	   sizeof(float) * bed->expCount);
    }
return(newBed);
}


int bedParseRgb(char *itemRgb)
/*      parse a string: "r,g,b" into three unsigned char values
        returned as 24 bit number, or -1 for failure */
{
char dupe[64];
int wordCount;
char *row[4];

strncpy(dupe, itemRgb, sizeof(dupe));
wordCount = chopString(dupe, ",", row, ArraySize(row));

if ((wordCount != 3) || (!isdigit(row[0][0]) ||
    !isdigit(row[1][0]) || !isdigit(row[2][0])))
        return (-1);

return ( ((atoi(row[0]) & 0xff) << 16) |
        ((atoi(row[1]) & 0xff) << 8) |
        (atoi(row[2]) & 0xff) );
}


void bedIntoRangeTree(struct bed *bed, struct rbTree *rangeTree)
/* Add all blocks in bed to range tree.  For beds without blocks,
 * add entire bed. */
{
int i;
if (bed->blockCount == 0)
    rangeTreeAdd(rangeTree, bed->chromStart, bed->chromEnd);
else
    {
    for (i=0; i < bed->blockCount; ++i)
	{
	int start = bed->chromStart + bed->chromStarts[i];
	int end = start + bed->blockSizes[i];
	rangeTreeAdd(rangeTree, start, end);
	}
    }
}

struct rbTree *bedToRangeTree(struct bed *bed)
/* Convert bed into a range tree. */
{
struct rbTree *rangeTree = rangeTreeNew();
bedIntoRangeTree(bed, rangeTree);
return rangeTree;
}

int bedRangeTreeOverlap(struct bed *bed, struct rbTree *rangeTree)
/* Return number of bases bed overlaps with rangeTree. */
{
int totalOverlap = 0;
if (bed->blockCount == 0)
    totalOverlap = rangeTreeOverlapSize(rangeTree, bed->chromStart, bed->chromEnd);
else
    {
    int i;
    for (i=0; i < bed->blockCount; ++i)
	{
	int start = bed->chromStart + bed->chromStarts[i];
	int end = start + bed->blockSizes[i];
	totalOverlap += rangeTreeOverlapSize(rangeTree, start, end);
	}
    }
return totalOverlap;
}

int bedSameStrandOverlap(struct bed *a, struct bed *b)
/* Return amount of block-level overlap on same strand between a and b */
{
/* Make sure on same strand, chromosome, and that overlap
 * at the non-block level. */
if (a->strand[0] != b->strand[0])
    return 0;
if (!sameString(a->chrom, b->chrom))
    return 0;
int outerOverlap = rangeIntersection(a->chromStart, a->chromEnd, 
	b->chromStart, b->chromEnd);
if (outerOverlap <= 0)
    return 0;

/* If both beds are non-blocked then we're pretty much done. */
if (a->blockCount == 0 && b->blockCount == 0)
    return outerOverlap;

/* Otherwise make up a range tree containing regions covered by a,
 * and figure out how much b overlaps it.. */
struct rbTree *rangeTree = bedToRangeTree(a);
int overlap = bedRangeTreeOverlap(b, rangeTree);

/* Clean up and return result. */
rangeTreeFree(&rangeTree);
return overlap;
}


struct bed *bedThickOnly(struct bed *in)
/* Return a bed that only has the thick part. (Which is usually the CDS). */
{
if (in->thickStart >= in->thickEnd)
     return NULL;
if (in->expCount != 0 || in->expIds != NULL || in->expScores != NULL)
   errAbort("Sorry, bedThickOnly only works on beds with up to 12 fields.");

/* Allocate output, and fill in simple fields. */
struct bed *out;
AllocVar(out);
out->chrom = cloneString(in->chrom);
out->chromStart = out->thickStart = in->thickStart;
out->chromEnd = out->thickEnd = in->thickEnd;
out->name = cloneString(in->name);
out->strand[0] = in->strand[0];
out->score = in->score;
out->itemRgb = in->itemRgb;

/* If need be fill in blocks. */
if (in->blockCount > 0)
   {
   /* Count up blocks inside CDS. */
   int i;
   int outBlockCount = 0;
   for (i=0; i<in->blockCount; ++i)
       {
       int start = in->chromStart + in->chromStarts[i];
       int end = start + in->blockSizes[i];
       if (start < in->thickStart) start = in->thickStart;
       if (end > in->thickEnd) end = in->thickEnd;
       if (start < end)
           outBlockCount += 1;
	}

    /* This trivieal case shouldn't happen, but just in case, we deal with it. */
    if (outBlockCount == 0)
        {
	freeMem(out);
	return NULL;
	}

    /* Allocate block arrays for output. */
    out->blockCount = outBlockCount;
    AllocArray(out->chromStarts, outBlockCount);
    AllocArray(out->blockSizes, outBlockCount);

    /* Scan through input one more time, copying to out. */
    int outBlockIx = 0;
    for (i=0; i<in->blockCount; ++i)
	{
	int start = in->chromStart + in->chromStarts[i];
	int end = start + in->blockSizes[i];
	if (start < in->thickStart) start = in->thickStart;
	if (end > in->thickEnd) end = in->thickEnd;
	if (start < end)
	    {
	    out->chromStarts[outBlockIx] = start - out->chromStart;
	    out->blockSizes[outBlockIx] = end - start;
	    outBlockIx += 1;
	    }
	}
    }
return out;
}


char *bedAsDef(int bedFieldCount, int totalFieldCount)
/* Return an autoSql definition for a bed of given number of fields. 
 * Normally totalFieldCount is equal to bedFieldCount.  If there are extra
 * fields they are just given the names field16, field17, etc and type string. */
{
if (bedFieldCount < 3 || bedFieldCount > 15)
    errAbort("bedFieldCount is %d, but must be between %d and %d in bedAsDef", bedFieldCount, 3, 15);
struct dyString *dy = dyStringNew(0);
dyStringAppend(dy, 
    "table bed\n"
    "\"Browser Extensible Data\"\n"
    "   (\n"
    "   string chrom;       \"Reference sequence chromosome or scaffold\"\n"
    "   uint   chromStart;  \"Start position in chromosome\"\n"
    "   uint   chromEnd;    \"End position in chromosome\"\n"
    );
if (bedFieldCount >= 4)
    dyStringAppend(dy, "   string name;        \"Name of item.\"\n");
if (bedFieldCount >= 5)
    dyStringAppend(dy, "   uint score;          \"Score (0-1000)\"\n");
if (bedFieldCount >= 6)
    dyStringAppend(dy, "   char[1] strand;     \"+ or - for strand\"\n");
if (bedFieldCount >= 7)
    dyStringAppend(dy, "   uint thickStart;   \"Start of where display should be thick (start codon)\"\n");
if (bedFieldCount >= 8)
    dyStringAppend(dy, "   uint thickEnd;     \"End of where display should be thick (stop codon)\"\n");
if (bedFieldCount >= 9)
    dyStringAppend(dy, "   uint reserved;     \"Used as itemRgb as of 2004-11-22\"\n");
if (bedFieldCount >= 10)
    dyStringAppend(dy, "   int blockCount;    \"Number of blocks\"\n");
if (bedFieldCount >= 11)
    dyStringAppend(dy, "   int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n");
if (bedFieldCount >= 12)
    dyStringAppend(dy, "   int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n");
if (bedFieldCount >= 13)
    dyStringAppend(dy, "   int expCount;	\"Experiment count\"\n");
if (bedFieldCount >= 14)
    dyStringAppend(dy, "   int[expCount] expIds;	\"Comma separated list of experiment ids. Always 0,1,2,3....\"\n");
if (bedFieldCount >= 15)
    dyStringAppend(dy, "   float[expCount] expScores; \"Comma separated list of experiment scores.\"\n");
int i;
for (i=bedFieldCount+1; i<=totalFieldCount; ++i)
    dyStringPrintf(dy, "lstring field%d;	\"Undocumented field\"\n", i+1);
dyStringAppend(dy, "   )\n");
return dyStringCannibalize(&dy);
}


void loadAndValidateBedExt(char *row[], int bedFieldCount, int fieldCount, struct lineFile *lf, struct bed * bed, struct asObject *as, boolean isCt,  boolean allow1bpOverlap)
/* Convert a row of strings to a bed and validate the contents.  Abort with message if invalid data. Optionally validate bedPlus via asObject.
 * If a customTrack, then some errors are tolerated. Possibly allow exons to overlap by one base. */
{
int count;
int *blockSizes = NULL;
int *chromStarts;

bed->chrom = row[0];  // note this value is not cloned for speed, callers may need to clone it.

// This check is usually redundant since the caller should be checking it against actual chromInfo names
// however hgLoadBed might not always have that info available.
if (strlen(bed->chrom) >= BB_MAX_CHROM_STRING)  // must leave room for 0 terminator
    lineFileAbort(lf, "chrom [%s] is too long (must not exceed %d characters)", bed->chrom, BB_MAX_CHROM_STRING - 1);
if (strlen(bed->chrom) < 1)
    lineFileAbort(lf, "chrom cannot be blank or empty");

lineFileAllInts(lf, row, 1, &bed->chromStart, FALSE, 4, "integer", FALSE);

lineFileAllInts(lf, row, 2, &bed->chromEnd, FALSE, 4, "integer", FALSE);

if (bed->chromEnd < bed->chromStart)
    lineFileAbort(lf, "chromStart after chromEnd (%u > %u)",
    	bed->chromStart, bed->chromEnd);
if (bedFieldCount > 3)
    {
    bed->name = row[3];
    if (strlen(bed->name) > 255)
	lineFileAbort(lf, "name [%s] is too long (must not exceed 255 characters)", bed->name);
    if (isCt)
	bed->name = cloneString(bed->name);
    }
if (bedFieldCount > 4)
    {
    lineFileAllInts(lf, row, 4, &bed->score, TRUE, 4, "integer", FALSE);
    if (!isCt && (bed->score < 0 || bed->score > 1000))
	    lineFileAbort(lf, "score (%d) must be between 0 and 1000", bed->score);
    }

if (bedFieldCount > 5)
    {
    if (!isCt && strlen(row[5]) > 1)
      lineFileAbort(lf, "Expecting + or - or . in strand, found [%s]",row[5]);
    bed->strand[0] = row[5][0];
    bed->strand[1] = 0;
    if (bed->strand[0] != '+' && bed->strand[0] != '-' && bed->strand[0] != '.')
      lineFileAbort(lf, "Expecting + or - or . in strand, found [%s]",row[5]);
    }
if (bedFieldCount > 6)
    lineFileAllInts(lf, row, 6, &bed->thickStart, FALSE, 4, "integer", FALSE);
else
    bed->thickStart = bed->chromStart;
if (bedFieldCount > 7)
    {
    lineFileAllInts(lf, row, 7, &bed->thickEnd, FALSE, 4, "integer", FALSE);
    if (bed->thickEnd < bed->thickStart)
     lineFileAbort(lf, "thickStart after thickEnd");
    if ((bed->thickStart != 0) &&
     ((bed->thickStart < bed->chromStart) ||
      (bed->thickStart > bed->chromEnd)))
     lineFileAbort(lf,
	 "thickStart out of range (chromStart to chromEnd, or 0 if no CDS)");
    if ((bed->thickEnd != 0) &&
     ((bed->thickEnd < bed->chromStart) ||
      (bed->thickEnd > bed->chromEnd)))
     lineFileAbort(lf,
	 "thickEnd out of range for %s:%u-%u, thick:%u-%u (chromStart to chromEnd, or 0 if no CDS)",
		   bed->name, bed->chromStart, bed->chromEnd,
		   bed->thickStart, bed->thickEnd);
    }
else
     bed->thickEnd = bed->chromEnd;

if (bedFieldCount > 8)
    {
    if (strchr(row[8],','))
	{
	unsigned char colors[4];
	char *saveColorString = cloneString(row[8]);
	int numColors = lineFileAllIntsArray(lf, row, 8, colors, sizeof colors, FALSE, 1, "integer", FALSE);
	if (numColors == 3)
	    {
	    bed->itemRgb = (((unsigned)colors[0]) << 2*8) | (((unsigned)colors[1]) << 1*8) | (unsigned)colors[2];
	    }
	else
	    lineFileAbort(lf, "Expecting color to consist of r,g,b values from 0 to 255. Got [%s]", saveColorString);
	freeMem(saveColorString);
	}
    else 
	{
	lineFileAllInts(lf, row, 8, &bed->itemRgb, FALSE, 4, "integer", FALSE);
	}
    }

int tempArraySize = 1;	// How big arrays are below
if (bedFieldCount > 9)
    {
    lineFileAllInts(lf, row, 9, &bed->blockCount, FALSE, 4, "integer", FALSE);
    if (!(bed->blockCount >= 1))
	lineFileAbort(lf, "Expecting blockCount (%d) to be 1 or more.", bed->blockCount);
    tempArraySize = bed->blockCount;
    }
int tempBlockSizes[tempArraySize];
int tempChromStarts[tempArraySize];
int tempExpIds[tempArraySize];
float tempExpScores[tempArraySize];
if (bedFieldCount > 10)
    {
    if (isCt)
	{
	AllocArray(bed->blockSizes,bed->blockCount+1); // having +1 allows us to detect incorrect size
        count = lineFileAllIntsArray(lf, row, 10, bed->blockSizes, bed->blockCount+1, TRUE, 4, "integer", TRUE);
	blockSizes = bed->blockSizes;
	}
    else
	{
        count = lineFileAllIntsArray(lf, row, 10, tempBlockSizes, tempArraySize, TRUE, 4, "integer", TRUE);
	blockSizes = tempBlockSizes;
	}
    if (count != bed->blockCount)
	lineFileAbort(lf, "Expecting %d elements in blockSizes list, found at least %d", bed->blockCount, count);
#ifdef NOTNOW
    int i;
    for (i=0; i < bed->blockCount;  i++)
	{
        if (!(blockSizes[i] > 0))
		lineFileAbort(lf, "BED blockSizes must be greater than 0, blockSize[%d] = %d", i, blockSizes[i]);
	}
#endif
    }
if (bedFieldCount > 11)
    {
    int i;
    if (isCt)
	{
	AllocArray(bed->chromStarts,bed->blockCount+1); // having +1 allows us to detect incorrect size
        count = lineFileAllIntsArray(lf, row, 11, bed->chromStarts, bed->blockCount+1, TRUE, 4, "integer", TRUE);
	chromStarts = bed->chromStarts;
	}
    else
	{
        count = lineFileAllIntsArray(lf, row, 11, tempChromStarts, tempArraySize, TRUE, 4, "integer", TRUE);
	chromStarts = tempChromStarts;
	}
    if (count != bed->blockCount)
	lineFileAbort(lf, "Expecting %d elements in chromStarts list, found at least %d", bed->blockCount, count);
    // tell the user if they appear to be using absolute starts rather than
    // relative... easy to forget!  Also check block order, coord ranges...
    if (chromStarts[0] != 0)
	lineFileAbort(lf,
	    "BED blocks must span chromStart to chromEnd.  "
	    "BED chromStarts[0] = %d, must be 0 so that (chromStart + "
	    "chromStarts[0]) equals chromStart.", chromStarts[0]);

    for (i=1; i < bed->blockCount;  i++)
	{

/*
printf("%d:%d %s %s s:%d c:%u cs:%u ce:%u csI:%d bsI:%d ls:%d le:%d<BR>\n", lineIx, i, bed->chrom, bed->name, bed->score, bed->blockCount, bed->chromStart, bed->chromEnd, bed->chromStarts[i], bed->blockSizes[i], lastStart, lastEnd);
*/
	// extra check to give user help for a common problem
	if (chromStarts[i]+bed->chromStart >= bed->chromEnd)
	    {
	    if (chromStarts[i] >= bed->chromStart)
		lineFileAbort(lf, "BED chromStarts offsets must be relative to chromStart, "
				  "not absolute.  Try subtracting chromStart from each offset "
				  "in chromStarts.");
	    else
		lineFileAbort(lf, "BED chromStarts[i]+chromStart must be less than chromEnd.");
	    }
	// chrom blocks must ascend without overlap
        int fudge = 0;
        if (allow1bpOverlap)
            fudge = -1;
        if (!(chromStarts[i] >= chromStarts[i-1] + blockSizes[i-1] + fudge))
		lineFileAbort(lf, "BED blocks must be in ascending order without overlap. Blocks %d and %d overlap.", i-1, i);
	}

    // last block-end must match chromEnd
    i = bed->blockCount-1;
    if ((bed->chromStart + chromStarts[i] + blockSizes[i]) != bed->chromEnd)
	{
	lineFileAbort(lf, "BED blocks must span chromStart to chromEnd.  (chromStart + "
			  "chromStarts[last] + blockSizes[last]) must equal chromEnd.");
	}
    }

if (bedFieldCount > 12)
    // get the microarray/colored-exon fields
    {
    lineFileAllInts(lf, row, 12, &bed->expCount, TRUE, 4, "integer", TRUE);
    if (!(bed->expCount >= 1))
	lineFileAbort(lf, "Expecting expCount (%d) to be 1 or more.", bed->expCount);
    if (isCt)
	{
	AllocArray(bed->expIds,bed->expCount+1); // having +1 allows us to detect incorrect size
        count = lineFileAllIntsArray(lf, row, 13, bed->expIds, bed->expCount+1, TRUE, 4, "integer", TRUE);
	}
    else
	{
        count = lineFileAllIntsArray(lf, row, 13, tempExpIds, tempArraySize, TRUE, 4, "integer", TRUE);
	}
    if (count != bed->expCount)
	lineFileAbort(lf, "expecting %d elements in expIds list (bed field 14)", bed->expCount);
    if (bedFieldCount == 15)
	{
	if (isCt)
	    {
	    sqlFloatDynamicArray(row[14], &bed->expScores, &count);
	    }
	else
	    {
	    count = sqlFloatArray(row[14], tempExpScores, tempArraySize);
	    }
	if (count != bed->expCount)
	    lineFileAbort(lf, "expecting %d elements in expScores list (bed field 15)", bed->expCount);
	}
    }

/* Check bedPlus fields are formatted right. */
/* This could form the basis of an .as-validator independent of BED. I suppose it could go in asParse.c */
if (as)
    {
    struct hash* linkHash = NULL;
    /* Validate as-fields */
    struct asColumn *asCol = NULL;
    asCol = as->columnList;
    int i;
    // Pre-scan ALL fields for linked fields
    for (i=0; i<fieldCount; ++i)
	{
	enum asTypes type = asCol->lowType->type;
	if (! (asCol->isList || asCol->isArray))
	    {
	    if (asTypesIsInt(type))
		{
		if (asCol->isSizeLink) // save the field value and index for later use in validating a list size.
		    {
		    int listSize = 0;  // big enough to hold the list count
		    lineFileAllInts(lf, row, i, &listSize, TRUE, 4, "integer", TRUE);
		    if (!linkHash)
			linkHash = newHash(4);
		    hashAddInt(linkHash, asCol->name, listSize);
		    }
		}
	    }
	asCol = asCol->next;
	}    
    /* Validate bed-plus fields */
    asCol = slElementFromIx(as->columnList, bedFieldCount);
    for (i=bedFieldCount; i<fieldCount; ++i)
	{
	enum asTypes type = asCol->lowType->type;
	if (! (asCol->isList || asCol->isArray))
	    {
	    if (asTypesIsInt(type))
		lineFileAllInts(lf, row, i, NULL, !asTypesIsUnsigned(type), asTypesIntSize(type), asTypesIntSizeDescription(type), FALSE);
	    else if (asTypesIsFloating(type))
		lineFileNeedDouble(lf, row, i);
	    else if (type == t_string)
		{
		if (strlen(row[i]) > 255)
		    lineFileAbort(lf, "expecting length (%llu) of string (%s) not to exceed 255 in field %s", (unsigned long long)strlen(row[i]), row[i], asCol->name);
		}
	    }
	else if (asCol->isList)
	    {
            if (asTypesIsFloating(type))
                {
                // assure count = #items in list; lightweight validation (better than none)
                int ix = asColumnFindIx(as->columnList, asCol->linkedSizeName);
                int count = sqlUnsigned(row[ix]);
		if (count < 0)
                    lineFileAbort(lf, 
                        "expecting nonnegative number in count field for %s list, found %d",
                                        asCol->name, asCol->fixedSize);
                int itemCount = countSeparatedItems(row[i], ',');
                if (count != itemCount)
                    lineFileAbort(lf, "expecting %d elements in %s list, found %d", 
                                        count, asCol->name, itemCount);
                }
	    else if (asTypesIsInt(type))
		{
		count = lineFileAllIntsArray(lf, row, i, NULL, countSeparatedItems(row[i], ','),
		    !asTypesIsUnsigned(type), asTypesIntSize(type), asTypesIntSizeDescription(type), FALSE);
		if (asCol->fixedSize > 0)
		    {
		    if (asCol->fixedSize != count)
			lineFileAbort(lf, "expecting %d elements in %s list, found %d", asCol->fixedSize, asCol->name, count);
		    }
		else
		    {
		    if (!linkHash)
			lineFileAbort(lf, "linked field %s was not found; it is required for determining listSize of %s"
			    , asCol->linkedSizeName, asCol->name);
		    int listSize = hashIntValDefault(linkHash, asCol->linkedSizeName, -1);
		    if (listSize == -1)
			lineFileAbort(lf, "linked field %s was not found; it is required for determining listSize of %s"
			    , asCol->linkedSizeName, asCol->name);
		    if (!(listSize >= 1))
			lineFileAbort(lf, "invalid list size %d for list %s must be 1 or greater, empty lists are not allowed", listSize, asCol->name);
		    if (!(listSize == count))
			lineFileAbort(lf, "expecting %d elements in %s list, found %d", listSize, asCol->name, count);
		    }
		}
	    }
	asCol = asCol->next;
	}
    hashFree(&linkHash);
    }

}


void bed3Free(struct bed3 **pBed)
/* Free up bed3 */
{
struct bed3 *bed = *pBed;
if (bed != NULL)
    {
    freeMem(bed->chrom);
    freez(pBed);
    }
}


void bed4Free(struct bed4 **pBed)
/* Free up bed4 */
{
struct bed4 *bed = *pBed;
if (bed != NULL)
    {
    freeMem(bed->chrom);
    freeMem(bed->name);
    freez(pBed);
    }
}
