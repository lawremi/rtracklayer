#ifndef BINRANGE_H
#define BINRANGE_H

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

#define BINRANGE_MAXEND_512M (512*1024*1024)
#define _binOffsetOldToExtended  4681


/*****  And now for some higher level stuff - useful for binning
 *****  things in memory. ******/

struct binElement
/* An element in a bin. */
    {
    struct binElement *next;
    int start, end;		/* 0 based, half open range */
    void *val;			/* Actual bin item. */
    };

int binElementCmpStart(const void *va, const void *vb);
/* Compare to sort based on start. */

struct binKeeper
/* This keeps things in bins in memory */
    {
    struct binKeeper *next;
    int minPos;		/* Minimum position to bin. */
    int maxPos;		/* Maximum position to bin. */
    int binCount;	/* Count of bins. */
    struct binElement **binLists; /* A list for each bin. */
    };

struct binKeeperCookie
/* used by binKeeperFirst/binKeeperNext in tracking location in traversing bins */
    {
    struct binKeeper *bk;       /* binKeeper we are associated with */
    int blIdx;                  /* current bin list index */
    struct binElement *nextBel; /* next binElement */
    };

struct binKeeper *binKeeperNew(int minPos, int maxPos);
/* Create new binKeeper that can cover range. */

void binKeeperAdd(struct binKeeper *bk, int start, int end, void *val);
/* Add item to binKeeper. */ 

struct binElement *binKeeperFind(struct binKeeper *bk, int start, int end);
/* Return a list of all items in binKeeper that intersect range.
 * Free this list with slFreeList. */

struct binElement *binKeeperFindSorted(struct binKeeper *bk, int start, int end);
/* Like binKeeperFind, but sort list on start coordinates. */

#endif /* BINRANGE_H */

