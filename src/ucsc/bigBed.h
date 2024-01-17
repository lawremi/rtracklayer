/* bigBed - interface to binary file with bed-style values (that is a bunch of
 * possibly overlapping regions.
 *
 * This shares a lot with the bigWig module. 
 *
 * Most of the functions here are concerned with reading bigBed files.  There's
 * two common things you want to do with a bigBed,  either stream through everything in it,
 * or just read the parts that overlap a interval within the genome.  The files
 * are optimized for interval queries, but streaming through them is not difficult either.
 * 
 * To query an interval:
 *    struct bbiFile *bbi = bigBedFileOpen(fileName);
 *    struct lm *lm = lmInit(0); // Memory pool to hold returned list
 *    struct bigBedInterval *list = bigBedIntervalQuery(bbi, chrom, start, end, 0, lm);
 *    struct bigBedInterval *el;
 *    for (el = list; el != NULL; el = el->next)
 *        // do something involving chrom, el->start, el->end
 *    lmCleanup(&lm);         // typically do this after each query
 *    bigBedFileClose(&bbi);  // typically only do this when finished all queries
 *
 * To stream through whole file
 *    struct bbiFile *bbi = bigBedFileOpen(fileName);
 *    struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
 *    for (chrom = chromList; chrom != NULL; chrom = chrom->next)
 *        {
 *        struct lm *lm = lmInit(0);
 *        struct bigBedInterval *list = bigBedIntervalQuery(bbi,chrom->name,0,chrom->size,0,lm);
 *        struct bigBedInterval *el;
 *        for (el = list; el != NULL; el = el->next)
 *            // do something involving chrom, el->start, el->end
 *        lmCleanup(&lm);
 *        }
 *    bigBedFileClose(&bbi);
 *
 * The processes for streaming through or doing interval queries on a bigWig file are very 
 * similar. */

#ifndef BIGBED_H
#define BIGBED_H

#ifndef BBIFILE
#include "bbiFile.h"
#endif

struct bigBedInterval
/* A partially parsed out bed record plus some extra fields.  Use this directly
 * or convert it to an array of characters with bigBedIntervalToRow. */
    {
    struct bigBedInterval *next;	/* Next in list. */
    bits32 start, end;		/* Range inside chromosome - half open zero based. */
    char *rest;			/* Rest of line. May be NULL*/
    bits32 chromId;             /* ID of chromosome.  */
    };

/*** Routines to open & close bigBed files, and to do chromosome range queries on them. ***/

struct bbiFile *bigBedFileOpen(char *fileName);
/* Open up big bed file.   Free this up with bbiFileClose. */

#define bigBedFileClose(a) bbiFileClose(a)

struct bigBedInterval *bigBedIntervalQuery(struct bbiFile *bbi, char *chrom,
	bits32 start, bits32 end, int maxItems, struct lm *lm);
/* Get data for interval.  Return list allocated out of lm.  Set maxItems to maximum
 * number of items to return, or to 0 for all items. */

int bigBedIntervalToRow(struct bigBedInterval *interval, char *chrom, char *startBuf, char *endBuf,
	char **row, int rowSize);
/* Convert bigBedInterval into an array of chars equivalent to what you'd get by
 * parsing the bed file. The startBuf and endBuf are used to hold the ascii representation of
 * start and end.  Note that the interval->rest string will have zeroes inserted as a side effect.
 * Returns number of fields in row.  */


/*** Some routines for accessing bigBed items via name. ***/


/** Routines to access other data from a bigBed file. */

char *bigBedAutoSqlText(struct bbiFile *bbi);
/* Get autoSql text if any associated with file.  Do a freeMem of this when done. */

#endif /* BIGBED_H */

