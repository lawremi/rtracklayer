/* rangeTree - This module is a way of keeping track of
 * non-overlapping ranges (half-open intervals). It is
 * based on the self-balancing rbTree code.  Use it in
 * place of a bitmap when the total number of ranges
 * is significantly smaller than the number of bits would
 * be. 
 * Beware the several static/global variables which can be
 * changed by various function calls. */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "limits.h"
#include "localmem.h"
#include "obscure.h"
#include "rbTree.h"
#include "rangeTree.h"


int rangeCmp(void *va, void *vb)
/* Return -1 if a before b,  0 if a and b overlap,
 * and 1 if a after b. */
{
struct range *a = va;
struct range *b = vb;
if (a->end <= b->start)
    return -1;
else if (b->end <= a->start)
    return 1;
else
    return 0;
}


static void *sumInt(void *a, void *b)
/* Local function used by rangeTreeAddValCount, which sums two ints a and b, 
 * referenced by void pointers, returning the result in a */
{
int *i = a, *j = b;
*i += *j;
return a;
}


struct range *rangeTreeAddVal(struct rbTree *tree, int start, int end, void *val, void *(*mergeVals)(void *existingVal, void *newVal) )
/* Add range to tree, merging with existing ranges if need be. 
 * If this is a new range, set the value to this val.
 * If there are existing items for this range, and if mergeVals function is not null, 
 * apply mergeVals to the existing values and this new val, storing the result as the val
 * for this range (see rangeTreeAddValCount() and rangeTreeAddValList() below for examples). */
{
struct range *r, *existing;
r = lmAlloc(tree->lm, sizeof(*r)); /* alloc new zeroed range */
r->start = start;
r->end = end;
r->val = val;
while ((existing = rbTreeRemove(tree, r)) != NULL)
    {
    r->start = min(r->start, existing->start);
    r->end = max(r->end, existing->end);
    if (mergeVals)
	r->val = mergeVals(existing->val, r->val);
    }
rbTreeAdd(tree, r);
return r;
}


struct range *rangeTreeAdd(struct rbTree *tree, int start, int end)
/* Add range to tree, merging with existing ranges if need be. */
{
    return rangeTreeAddVal(tree, start, end, NULL, NULL);
}


void rangeTreeAddToCoverageDepth(struct rbTree *tree, int start, int end)
/* Add area from start to end to a tree that is being built up to store the
 * depth of coverage.  Recover coverage back out by looking at ptToInt(range->val)
 * on tree elements. */
{
struct range q;
q.start = start;
q.end = end;

struct range *r, *existing = rbTreeFind(tree, &q);
if (existing == NULL)
    {
    lmAllocVar(tree->lm, r);
    r->start = start;
    r->end = end;
    r->val = intToPt(1);
    rbTreeAdd(tree, r);
    }
else
    {
    if (existing->start <= start && existing->end >= end)
    /* The existing one completely encompasses us */
        {
	/* Make a new section for the bit before start. */
	if (existing->start < start)
	    {
	    lmAllocVar(tree->lm, r);
	    r->start = existing->start;
	    r->end = start;
	    r->val = existing->val;
	    existing->start = start;
	    rbTreeAdd(tree, r);
	    }
	/* Make a new section for the bit after end. */
	if (existing->end > end)
	    {
	    lmAllocVar(tree->lm, r);
	    r->start = end;
	    r->end = existing->end;
	    r->val = existing->val;
	    existing->end = end;
	    rbTreeAdd(tree, r);
	    }
	/* Increment existing section in overlapping area. */
        existing->val = (char *)(existing->val) + 1;
	}
    else
    /* In general case fetch list of regions that overlap us. 
       Remaining cases to handle are: 
	     r >> e     rrrrrrrrrrrrrrrrrrrr
			     eeeeeeeeee

	     e < r           rrrrrrrrrrrrrrr
			eeeeeeeeeeee

	     r < e      rrrrrrrrrrrr
			     eeeeeeeeeeeee
     */
        {
	struct range *existingList = rangeTreeAllOverlapping(tree, start, end);

#ifdef DEBUG
	/* Make sure that list is really sorted for debugging... */
	int lastStart = existingList->start;
	for (r = existingList; r != NULL; r = r->next)
	    {
	    int start = r->start;
	    if (start < lastStart)
	        internalErr();
	    }
#endif /* DEBUG */

	int s = start, e = end;
	for (existing = existingList; existing != NULL; existing = existing->next)
	    {
	    /* Deal with start of new range that comes before existing */
	    if (s < existing->start)
	        {
		lmAllocVar(tree->lm, r);
		r->start = s;
		r->end = existing->start;
		r->val = intToPt(1);
		s = existing->start;
		rbTreeAdd(tree, r);
		}
	    else if (s > existing->start)
	        {
		lmAllocVar(tree->lm, r);
		r->start = existing->start;
		r->end = s;
		r->val = existing->val;
		existing->start = s;
		rbTreeAdd(tree, r);
		}
	    existing->val = (char *)(existing->val) + 1;
	    s = existing->end;
	    }
	if (s < e)
	/* Deal with end of new range that doesn't overlap with anything. */
	    {
	    lmAllocVar(tree->lm, r);
	    r->start = s;
	    r->end = e;
	    r->val = intToPt(1);
	    rbTreeAdd(tree, r);
	    }
	}
    }

}

static struct range *rangeList;

static void rangeListAdd(void *v)
/* Callback to add item to range list. */
{
struct range *r = v;
slAddHead(&rangeList, r);
}

struct range *rangeTreeList(struct rbTree *tree)
/* Return list of all ranges in tree in order.  Not thread safe. 
 * No need to free this when done, memory is local to tree. */
{
rangeList = NULL;
rbTreeTraverse(tree, rangeListAdd);
slReverse(&rangeList);
return rangeList;
}

struct range *rangeTreeAllOverlapping(struct rbTree *tree, int start, int end)
/* Return list of all items in range tree that overlap interval start-end.
 * Do not free this list, it is owned by tree.  However it is only good until
 * next call to rangeTreeFindInRange or rangeTreeList. Not thread safe. */
{
struct range tempR;
tempR.start = start;
tempR.end = end;
rangeList = NULL;
rbTreeTraverseRange(tree, &tempR, &tempR, rangeListAdd);
slReverse(&rangeList);
return rangeList;
}


/* A couple of variables used to calculate total overlap. */
static int totalOverlap;
static int overlapStart, overlapEnd;

static void addOverlap(void *v)
/* Callback to add item to range list. */
{
struct range *r = v;
totalOverlap += positiveRangeIntersection(r->start, r->end, 
	overlapStart, overlapEnd);
}

int rangeTreeOverlapSize(struct rbTree *tree, int start, int end)
/* Return the total size of intersection between interval
 * from start to end, and items in range tree. Sadly not
 * thread-safe. 
 * On 32 bit machines be careful not to overflow
 * range of start, end or total size return value. */
{
struct range tempR;
tempR.start = overlapStart = start;
tempR.end = overlapEnd = end;
totalOverlap = 0;
rbTreeTraverseRange(tree, &tempR, &tempR, addOverlap);
return totalOverlap;
}

void rangeTreeSumRangeCallback(void *item, void *context)
/* This is a callback for rbTreeTraverse with context.  It just adds up
 * end-start */
{
struct range *range = item;
long long *pSum = context;
*pSum += range->end - range->start;
}


struct rbTree *rangeTreeNew()
/* Create a new, empty, rangeTree. */
{
return rbTreeNew(rangeCmp);
}
