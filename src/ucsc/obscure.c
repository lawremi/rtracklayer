/* Obscure stuff that is handy every now and again. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include <unistd.h>
#include "portable.h"
#include "localmem.h"
#include "hash.h"
#include "obscure.h"
#include "linefile.h"

static int _dotForUserMod = 100; /* How often does dotForUser() output a dot. */

int digitsBaseTwo(unsigned long x)
/* Return base two # of digits. */
{
int digits = 0;
while (x)
    {
    digits += 1;
    x >>= 1;
    }
return digits;
}

void writeGulp(char *file, char *buf, int size)
/* Write out a bunch of memory. */
{
FILE *f = mustOpen(file, "w");
mustWrite(f, buf, size);
carefulClose(&f);
}

void readInGulp(char *fileName, char **retBuf, size_t *retSize)
/* Read whole file in one big gulp. */
{
size_t size = (size_t)fileSize(fileName);
char *buf;
FILE *f = mustOpen(fileName, "rb");
*retBuf = buf = needLargeMem(size+1);
mustRead(f, buf, size);
buf[size] = 0;      /* Just in case it needs zero termination. */
fclose(f);
if (retSize != NULL)
    *retSize = size;
}

void *intToPt(int i)
/* Convert integer to pointer. Use when really want to store an
 * int in a pointer field. */
{
char *pt = NULL;
return pt+i;
}

int ptToInt(void *pt)
/* Convert pointer to integer.  Use when really want to store a
 * pointer in an int. */
{
char *a = NULL, *b = pt;
return b - a;
}

boolean parseQuotedString( char *in, char *out, char **retNext)
/* Read quoted string from in (which should begin with first quote).
 * Write unquoted string to out, which may be the same as in.
 * Return pointer to character past end of string in *retNext. 
 * Return FALSE if can't find end. */
{
char c, *s = in;
int quoteChar = *s++;
boolean escaped = FALSE;

for (;;)
   {
   c = *s++;
   if (c == 0)
       {
       warn("Unmatched %c", quoteChar);
       return FALSE;
       }
   if (escaped)
       {
       if (c == '\\' || c == quoteChar)
          *out++ = c;
       else
          {
	  *out++ = '\\';
	  *out++ = c;
	  }
       escaped = FALSE;
       }
   else
       {
       if (c == '\\')
           escaped = TRUE;
       else if (c == quoteChar)
           break;
       else
           *out++ = c;
       }
   }
*out = 0;
if (retNext != NULL)
    *retNext = s;
return TRUE;
}

void escCopy(char *in, char *out, char toEscape, char escape)
/* Copy in to out, escaping as needed.  Out better be big enough. 
 * (Worst case is strlen(in)*2 + 1.) */
{
char c;
for (;;)
    {
    c = *in++;
    if (c == toEscape)
        *out++ = escape;
    *out++ = c;
    if (c == 0)
        break;
    }
}

struct hash *hashThisEqThatLine(char *line, int lineIx, boolean firstStartsWithLetter)
/* Return a symbol table from a line of form:
 *   1-this1=val1 2-this='quoted val2' var3="another val" 
 * If firstStartsWithLetter is true, then the left side of the equals must start with
 * a letter. */
{
char *dupe = cloneString(line);
char *s = dupe, c;
char *var, *val;
struct hash *hash = newHash(8);

for (;;)
    {
    if ((var = skipLeadingSpaces(s)) == NULL)
        break;

    if ((c = *var) == 0)
        break;
    if (firstStartsWithLetter && !isalpha(c))
	errAbort("line %d of custom input: variable needs to start with letter '%s'", lineIx, var);
    val = strchr(var, '=');
    if (val == NULL)
        {
        errAbort("line %d of var %s in custom input: %s \n missing = in var/val pair", lineIx, var, line);
        }
    *val++ = 0;
    c = *val;
    if (c == '\'' || c == '"')
        {
	if (!parseQuotedString(val, val, &s))
	    errAbort("line %d of input: missing closing %c", lineIx, c);
	}
    else
	{
	s = skipToSpaces(val);
	if (s != NULL) *s++ = 0;
	}
    hashAdd(hash, var, cloneString(val));
    }
freez(&dupe);
return hash;
}


struct slName *charSepToSlNames(char *string, char c)
/* Convert character-separated list of items to slName list. 
 * Note that the last occurence of c is optional.  (That
 * is for a comma-separated list a,b,c and a,b,c, are
 * equivalent. */
{
struct slName *list = NULL, *el;
char *s, *e;

s = string;
while (s != NULL && s[0] != 0)
    {
    e = strchr(s, c);
    if (e == NULL)
        {
	el = slNameNew(s);
	slAddHead(&list, el);
	break;
	}
    else
        {
	el = slNameNewN(s, e - s);
	slAddHead(&list, el);
	s = e+1;
	}
    }
slReverse(&list);
return list;
}


void sprintLongWithCommas(char *s, long long l)
/* Print out a long number with commas a thousands, millions, etc. */
{
long long trillions, billions, millions, thousands;
if (l >= 1000000000000LL)
    {
    trillions = l/1000000000000LL;
    l -= trillions * 1000000000000LL;
    billions = l/1000000000;
    l -= billions * 1000000000;
    millions = l/1000000;
    l -= millions * 1000000;
    thousands = l/1000;
    l -= thousands * 1000;
    sprintf(s, "%lld,%03lld,%03lld,%03lld,%03lld", trillions, billions, millions, thousands, l);
    }
else if (l >= 1000000000)
    {
    billions = l/1000000000;
    l -= billions * 1000000000;
    millions = l/1000000;
    l -= millions * 1000000;
    thousands = l/1000;
    l -= thousands * 1000;
    sprintf(s, "%lld,%03lld,%03lld,%03lld", billions, millions, thousands, l);
    }
else if (l >= 1000000)
    {
    millions = l/1000000;
    l -= millions * (long long)1000000;
    thousands = l/1000;
    l -= thousands * 1000;
    sprintf(s, "%lld,%03lld,%03lld", millions, thousands, l);
    }
else if (l >= 1000)
    {
    thousands = l/1000;
    l -= thousands * 1000;
    sprintf(s, "%lld,%03lld", thousands, l);
    }
else
    sprintf(s, "%lld", l);
}

void shuffleArrayOfPointers(void *pointerArray, int arraySize)
/* Shuffle array of pointers of given size given number of times. */
{
void **array = pointerArray, *pt;
int i, randIx;

/* Randomly permute an array using the method from Cormen, et al */
for (i=0; i<arraySize; ++i)
    {
    randIx = i + (rand() % (arraySize - i));
    pt = array[i];
    array[i] = array[randIx];
    array[randIx] = pt;
    }
}

void shuffleList(void *pList)
/* Randomize order of slList.  Usage:
 *     randomizeList(&list)
 * where list is a pointer to a structure that
 * begins with a next field. */
{
struct slList **pL = (struct slList **)pList;
struct slList *list = *pL;
int count;
count = slCount(list);
if (count > 1)
    {
    struct slList *el;
    struct slList **array;
    int i;
    array = needLargeMem(count * sizeof(*array));
    for (el = list, i=0; el != NULL; el = el->next, i++)
        array[i] = el;
    for (i=0; i<4; ++i)
        shuffleArrayOfPointers(array, count);
    list = NULL;
    for (i=0; i<count; ++i)
        {
        array[i]->next = list;
        list = array[i];
        }
    freeMem(array);
    slReverse(&list);
    *pL = list;       
    }
}

void *slListRandomReduce(void *list, double reduceRatio)
/* Reduce list to approximately reduceRatio times original size. Destroys original list. */
{
if (reduceRatio >= 1.0)
    return list;
int threshold = RAND_MAX * reduceRatio;
struct slList *newList = NULL, *next, *el;
for (el = list; el != NULL; el = next)
    {
    next = el->next;
    if (rand() <= threshold)
        {
	slAddHead(&newList, el);
	}
    }
return newList;
}
