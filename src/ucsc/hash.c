/* Hash.c - implements hashing.  See hash.h for usage comments.
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "localmem.h"
#include "hash.h"
#include "obscure.h"
#include "dystring.h"


/*
 * Hash a string key.  This code is taken from Tcl interpreter. I was borrowed
 * after discovering a lot of collisions and poor utilization of the table
 * when hashing accessions.
 *
 * This function was compared to Bob Jenkins' lookup2 hash function and
 * (http://burtleburtle.net/bob/hash/) and Paul Hsieh's SuperFast
 * hash function (http://www.azillionmonkeys.com/qed/hash.html).
 * Both of those functions provided better utilization of the table,
 * but were also more expensive, so the Tcl function was used.
 * If hashing of binary keys is implemented, SuperFast hash should
 * be considered.
 *
 * for an explanation of this function, see HashStringKey() in the
 * Tcl source file, generic/tclHash.c, available from
 * http://tcl.sourceforge.net/.
 *
 * The Tcl code is:
 * Copyright (c) 1991-1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
 *
 * See the file "license.terms" (in the Tcl distribution) for complete
 * license (which is a BSD-style license).
 *
 * Since hashCrc() is in use elsewhere, 
 * a new function hashString() was created for use in hash table.
 * -- markd
 */
bits32 hashString(char *string)
/* Compute a hash value of a string. */
{
char *keyStr = string;
unsigned int result = 0;
int c;

while ((c = *keyStr++) != '\0')
    {
    result += (result<<3) + c;
    }
return result;
}

struct hashEl *hashLookup(struct hash *hash, char *name)
/* Looks for name in hash table. Returns associated element,
 * if found, or NULL if not.  If there are multiple entries
 * for name, the last one added is returned (LIFO behavior).
 */
{
struct hashEl *el = hash->table[hashString(name)&hash->mask];
while (el != NULL)
    {
    if (strcmp(el->name, name) == 0)
        break;
    el = el->next;
    }
return el;
}

struct hashEl *hashLookupUpperCase(struct hash *hash, char *name)
/* Lookup upper cased name in hash. (Assumes all elements of hash
 * are themselves already in upper case.) */
{
char s[256];
safef(s, sizeof(s), "%s", name);
touppers(s);
return hashLookup(hash, s);
}

struct hashEl *hashAddN(struct hash *hash, char *name, int nameSize, void *val)
/* Add name of given size to hash (no need to be zero terminated) */
{
struct hashEl *el;
if (hash->lm) 
    el = lmAlloc(hash->lm, sizeof(*el));
else
    AllocVar(el);
el->hashVal = hashString(name);
int hashVal = el->hashVal & hash->mask;
if (hash->lm)
    {
    el->name = lmAlloc(hash->lm, nameSize+1);
    memcpy(el->name, name, nameSize);
    }
else
    el->name = cloneStringZ(name, nameSize);
el->val = val;
el->next = hash->table[hashVal];
hash->table[hashVal] = el;
hash->elCount += 1;
if (hash->autoExpand && hash->elCount > (int)(hash->size * hash->expansionFactor))
    {
    /* double the size */
    hashResize(hash, digitsBaseTwo(hash->size));
    }
return el;
}

struct hashEl *hashAdd(struct hash *hash, char *name, void *val)
/* Add new element to hash table.  If an item with name, already exists, a new
 * item is added in a LIFO manner.  The last item added for a given name is
 * the one returned by the hashLookup functions.  hashLookupNext must be used
 * to find the preceding entries for a name.
 */
{
return hashAddN(hash, name, strlen(name), val);
}

void freeHashEl(struct hashEl *hel)
/* Free hash element. Use only on non-local memory version. */
{
freeMem(hel->name);
freeMem(hel);
}

void *hashRemove(struct hash *hash, char *name)
/* Remove item of the given name from hash table. 
 * Returns value of removed item, or NULL if not in the table.
 * If their are multiple entries for name, the last one added
 * is removed (LIFO behavior).
 */
{
struct hashEl *hel;
void *ret;
struct hashEl **pBucket = &hash->table[hashString(name)&hash->mask];
for (hel = *pBucket; hel != NULL; hel = hel->next)
    if (sameString(hel->name, name))
        break;
if (hel == NULL)
    return NULL;
ret = hel->val;
if (slRemoveEl(pBucket, hel))
    {
    hash->elCount -= 1;
    if (!hash->lm)
	freeHashEl(hel);
    }
return ret;
}

struct hashEl *hashAddUnique(struct hash *hash, char *name, void *val)
/* Add new element to hash table. Squawk and die if not unique */
{
if (hashLookup(hash, name) != NULL)
    errAbort("%s duplicated, aborting", name);
return hashAdd(hash, name, val);
}

struct hashEl *hashAddSaveName(struct hash *hash, char *name, void *val, char **saveName)
/* Add new element to hash table.  Save the name of the element, which is now
 * allocated in the hash table, to *saveName.  A typical usage would be:
 *    AllocVar(el);
 *    hashAddSaveName(hash, name, el, &el->name);
 */
{
struct hashEl *hel = hashAdd(hash, name, val);
*saveName = hel->name;
return hel;
}

int hashIntVal(struct hash *hash, char *name)
/* Return integer value associated with name in a simple 
 * hash of ints. */
{
void *val = hashMustFindVal(hash, name);
return ptToInt(val);
}

int hashIntValDefault(struct hash *hash, char *name, int defaultInt)
/* Return integer value associated with name in a simple 
 * hash of ints or defaultInt if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if(hel == NULL)
    return defaultInt;
return ptToInt(hel->val);
}

void *hashMustFindVal(struct hash *hash, char *name)
/* Lookup name in hash and return val.  Abort if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
    errAbort("hashMustFindVal: '%s' not found", name);
return hel->val;
}

void *hashFindVal(struct hash *hash, char *name)
/* Look up name in hash and return val or NULL if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
    return NULL;
return hel->val;
}

void *hashFindValUpperCase(struct hash *hash, char *name)
/* Lookup upper cased name in hash and return val or return NULL if not found.
 * (Assumes all elements of hash are themselves already in upper case.) */
{
struct hashEl *hel = hashLookupUpperCase(hash, name);
if (hel == NULL)
    return NULL;
return hel->val;
}

struct hashEl *hashAddInt(struct hash *hash, char *name, int val)
/* Store integer value in hash */
{
char *pt = NULL;
return hashAdd(hash, name, pt + val);
}

struct hash *newHashExt(int powerOfTwoSize, boolean useLocalMem)
/* Returns new hash table. Uses local memory optionally. */
{
struct hash *hash = needMem(sizeof(*hash));
int memBlockPower = 16;
if (powerOfTwoSize == 0)
    powerOfTwoSize = 12;
assert(powerOfTwoSize <= hashMaxSize && powerOfTwoSize > 0);
hash->powerOfTwoSize = powerOfTwoSize;
hash->size = (1<<powerOfTwoSize);
/* Make size of memory block for allocator vary between
 * 256 bytes and 64k depending on size of table. */
if (powerOfTwoSize < 8)
    memBlockPower = 8;
else if (powerOfTwoSize < 16)
    memBlockPower = powerOfTwoSize;
if (useLocalMem) 
    hash->lm = lmInit(1<<memBlockPower);
hash->mask = hash->size-1;
AllocArray(hash->table, hash->size);
hash->autoExpand = TRUE;
hash->expansionFactor = defaultExpansionFactor;   /* Expand when elCount > size*expansionFactor */
return hash;
}

void hashResize(struct hash *hash, int powerOfTwoSize)
/* Resize the hash to a new size */
{
int oldHashSize = hash->size;
struct hashEl **oldTable = hash->table;

if (powerOfTwoSize == 0)
    powerOfTwoSize = 12;
assert(powerOfTwoSize <= hashMaxSize && powerOfTwoSize > 0);
hash->powerOfTwoSize = powerOfTwoSize;
hash->size = (1<<powerOfTwoSize);
hash->mask = hash->size-1;

AllocArray(hash->table, hash->size);

int i;
struct hashEl *hel, *next;
for (i=0; i<oldHashSize; ++i)
    {
    for (hel = oldTable[i]; hel != NULL; hel = next)
	{
	next = hel->next;
	int hashVal = hel->hashVal & hash->mask;
	hel->next = hash->table[hashVal];
	hash->table[hashVal] = hel;
	}
    }
/* restore original list order */
for (i=0; i<hash->size; ++i)
    {
    struct hashEl *hel = hash->table[i];
    if (hel != NULL && hel->next != NULL)	    
	slReverse(&hash->table[i]);
    }
freeMem(oldTable);
hash->numResizes++;
}


void hashTraverseVals(struct hash *hash, void (*func)(void *val))
/* Apply func to every element of hash with hashEl->val as parameter. */
{
int i;
struct hashEl *hel;
for (i=0; i<hash->size; ++i)
    {
    for (hel = hash->table[i]; hel != NULL; hel = hel->next)
	func(hel->val);
    }
}

int hashElCmp(const void *va, const void *vb)
/* Compare two hashEl by name. */
{
const struct hashEl *a = *((struct hashEl **)va);
const struct hashEl *b = *((struct hashEl **)vb);
return strcmp(a->name, b->name);
}

struct hashEl *hashElListHash(struct hash *hash)
/* Return a list of all elements of hash.   Free return with hashElFreeList. */
{
int i;
struct hashEl *hel, *dupe, *list = NULL;
for (i=0; i<hash->size; ++i)
    {
    for (hel = hash->table[i]; hel != NULL; hel = hel->next)
	{
	dupe = CloneVar(hel);
	slAddHead(&list, dupe);
	}
    }
return list;
}


void hashElFreeList(struct hashEl **pList)
/* Free hash el list returned from hashListAll.  (Don't use
 * this internally. */
{
slFreeList(pList);
}

struct hashCookie hashFirst(struct hash *hash)
/* Return an object to use by hashNext() to traverse the hash table.
 * The first call to hashNext will return the first entry in the table. */
{
struct hashCookie cookie;
cookie.hash = hash;
cookie.idx = 0;
cookie.nextEl = NULL;

/* find first entry */
for (cookie.idx = 0;
     (cookie.idx < hash->size) && (hash->table[cookie.idx] == NULL);
     cookie.idx++)
    continue;  /* empty body */
if (cookie.idx < hash->size)
    cookie.nextEl = hash->table[cookie.idx];
return cookie;
}

struct hashEl* hashNext(struct hashCookie *cookie)
/* Return the next entry in the hash table, or NULL if no more. Do not modify
 * hash table while this is being used. */
{
/* NOTE: if hashRemove were coded to track the previous entry during the
 * search and then use it to do the remove, it would be possible to
 * remove the entry returned by this method */
struct hashEl *retEl = cookie->nextEl;
if (retEl == NULL)
    return NULL;  /* no more */

/* find next entry */
cookie->nextEl = retEl->next;
if (cookie->nextEl == NULL)
    {
    for (cookie->idx++; (cookie->idx < cookie->hash->size)
             && (cookie->hash->table[cookie->idx] == NULL); cookie->idx++)
        continue;  /* empty body */
    if (cookie->idx < cookie->hash->size)
        cookie->nextEl = cookie->hash->table[cookie->idx];
    }
return retEl;
}

void freeHash(struct hash **pHash)
/* Free up hash table. */
{
struct hash *hash = *pHash;
if (hash == NULL)
    return;
if (hash->lm)
    lmCleanup(&hash->lm);
else
    {
    int i;
    struct hashEl *hel, *next;
    for (i=0; i<hash->size; ++i)
	{
	for (hel = hash->table[i]; hel != NULL; hel = next)
	    {
	    next = hel->next;
	    freeHashEl(hel);
	    }
	}
    }
freeMem(hash->table);
freez(pHash);
}


static int bucketLen(struct hashEl *hel)
/* determine how many elements are in a hash bucket */
{
int nel = 0;
for (; hel != NULL; hel = hel->next)
    nel++;
return nel;
}

int hashNumEntries(struct hash *hash)
/* count the number of entries in a hash */
{
int n = 0, i;
for (i=0; i<hash->size; ++i)
    n += bucketLen(hash->table[i]);
return n;
}
