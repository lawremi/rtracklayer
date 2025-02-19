/* memalloc.c - Routines to allocate and deallocate dynamic memory. 
 * This lets you have a stack of memory handlers.  The default
 * memory handler is a thin shell around malloc/free.  You can
 * substitute routines that do more integrety checking with
 * pushCarefulMem(), or routines of your own devising with
 * pushMemHandler(). 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include <pthread.h>
#include "common.h"
#include "obscure.h"
#include "memalloc.h"
#include "dlist.h"


static void *defaultAlloc(size_t size)
/* Default allocator. */
{
return malloc(size);
}

static void defaultFree(void *vpt)
/* Default deallocator. */
{
free(vpt);
}

static void *defaultRealloc(void *vpt, size_t size)
/* Default deallocator. */
{
return realloc(vpt, size);
}

static struct memHandler defaultMemHandler = 
/* Default memory handler. */
    {
    NULL,
    defaultAlloc,
    defaultFree,
    defaultRealloc,
    };

static struct memHandler *mhStack = &defaultMemHandler;

struct memHandler *pushMemHandler(struct memHandler *newHandler)
/* Use newHandler for memory requests until matching popMemHandler.
 * Returns previous top of memory handler stack. */
{
struct memHandler *oldHandler = mhStack;
slAddHead(&mhStack, newHandler);
return oldHandler;
}


struct memHandler *popMemHandler()
/* Removes top element from memHandler stack and returns it. */
{
struct memHandler *oldHandler = mhStack;
if (mhStack == &defaultMemHandler)
    errAbort("Too many popMemHandlers()");
mhStack = mhStack->next;
return oldHandler;
}


/* 128*8*1024*1024 == 1073741824 == 2^30 on 32 bit machines,size_t == 4 bytes*/
/* on 64 bit machines, size_t = 8 bytes, 2^30 * 2 * 2 * 2 * 2 = 2^34 == 16 Gb */
static size_t maxAlloc = (size_t)128*8*1024*1024*(sizeof(size_t)/4)*(sizeof(size_t)/4)*(sizeof(size_t)/4*(sizeof(size_t)/4));

void *needLargeMem(size_t size)
/* This calls abort if the memory allocation fails. The memory is
 * not initialized to zero. */
{
void *pt;
if (size == 0 || size >= maxAlloc)
    errAbort("needLargeMem: trying to allocate %llu bytes (limit: %llu)",
         (unsigned long long)size, (unsigned long long)maxAlloc);
if ((pt = mhStack->alloc(size)) == NULL)
    errAbort("needLargeMem: Out of memory - request size %llu bytes, errno: %d\n",
             (unsigned long long)size, errno);
return pt;
}

void *needLargeZeroedMem(size_t size)
/* Request a large block of memory and zero it. */
{
void *v;
v = needLargeMem(size);
memset(v, 0, size);
return v;
}

void *needLargeMemResize(void* vp, size_t size)
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL,
 * a new memory block is allocated.  Memory not initted. */
{
void *pt;
if (size == 0 || size >= maxAlloc)
    errAbort("needLargeMemResize: trying to allocate %llu bytes (limit: %llu)",
         (unsigned long long)size, (unsigned long long)maxAlloc);
if ((pt = mhStack->realloc(vp, size)) == NULL)
    errAbort("needLargeMemResize: Out of memory - request size %llu bytes, errno: %d\n",
             (unsigned long long)size, errno);
return pt;
}

void *needLargeZeroedMemResize(void* vp, size_t oldSize, size_t newSize)
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL, a
 * new memory block is allocated.  If block is grown, new memory is zeroed. */
{
void *v = needLargeMemResize(vp, newSize);
if (newSize > oldSize)
    memset(((char*)v)+oldSize, 0, newSize-oldSize);
return v;
}

void *needHugeMem(size_t size)
/* No checking on size.  Memory not initted. */
{
void *pt;
if (size == 0)
    errAbort("needHugeMem: trying to allocate 0 bytes");
if ((pt = mhStack->alloc(size)) == NULL)
    errAbort("needHugeMem: Out of huge memory - request size %llu bytes, errno: %d\n",
             (unsigned long long)size, errno);
return pt;
}


void *needHugeMemResize(void* vp, size_t size)
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL,
 * a new memory block is allocated.  No checking on size.  Memory not
 * initted. */
{
void *pt;
if ((pt = mhStack->realloc(vp, size)) == NULL)
    errAbort("needHugeMemResize: Out of memory - request resize %llu bytes, errno: %d\n",
	(unsigned long long)size, errno);
return pt;
}


#define NEEDMEM_LIMIT 500000000

void *needMem(size_t size)
/* Need mem calls abort if the memory allocation fails. The memory
 * is initialized to zero. */
{
void *pt;
if (size == 0 || size > NEEDMEM_LIMIT)
    errAbort("needMem: trying to allocate %llu bytes (limit: %llu)",
         (unsigned long long)size, (unsigned long long)NEEDMEM_LIMIT);
if ((pt = mhStack->alloc(size)) == NULL)
    errAbort("needMem: Out of memory - request size %llu bytes, errno: %d\n",
             (unsigned long long)size, errno);
memset(pt, 0, size);
return pt;
}

void *needMoreMem(void *old, size_t oldSize, size_t newSize)
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL, a
 * new memory block is allocated.  No checking on size.  If block is grown,
 * new memory is zeroed. */
{
return needLargeZeroedMemResize(old, oldSize, newSize);
}

void freeMem(void *pt)
/* Free memory will check for null before freeing. */
{
if (pt != NULL)
    mhStack->free(pt);
}

void freez(void *vpt)
/* Pass address of pointer.  Will free pointer and set it 
 * to NULL. */
{
void **ppt = (void **)vpt;
void *pt = *ppt;
*ppt = NULL;
freeMem(pt);
}
