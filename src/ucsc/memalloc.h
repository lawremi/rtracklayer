/* Let the user redirect where memory allocation/deallocation
 * happens.  'careful' routines help debug scrambled heaps. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef MEMALLOC_H
#define MEMALLOC_H

struct memHandler
    {
    struct memHandler *next;
    void * (*alloc)(size_t size);
    void (*free)(void *vpt);
    void * (*realloc)(void* vpt, size_t size);
    };

struct memHandler *pushMemHandler(struct memHandler *newHandler);
/* Use newHandler for memory requests until matching popMemHandler.
 * Returns previous top of memory handler stack. */

struct memHandler *popMemHandler();
/* Removes top element from memHandler stack and returns it. */

#endif /* MEMALLOC_H */

