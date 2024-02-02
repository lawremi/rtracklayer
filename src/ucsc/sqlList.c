/* Stuff for processing comma separated lists - a little long so
 * in a separate module from jksql.c though interface is still
 * in jksql.c. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

/* The various static routines sql<Type>StaticArray are NOT thread-safe. */

#include "common.h"
#include "sqlNum.h"
#include "sqlList.h"
#include "dystring.h"
#include "hash.h"

int sqlDoubleArray(char *s, double *array, int maxArraySize)
/* Convert comma separated list of floating point numbers to an array.  
 * Pass in array and max size of array. */
{
unsigned count = 0;
for (;;)
    {
    char *e;
    if (s == NULL || s[0] == 0 || count == maxArraySize)
	break;
    e = strchr(s, ',');
    if (e != NULL)
	*e++ = 0;
    array[count++] = atof(s);
    s = e;
    }
return count;
}

int sqlFloatArray(char *s, float *array, int maxArraySize)
/* Convert comma separated list of floating point numbers to an array.  
 * Pass in array and max size of array. */
{
unsigned count = 0;
for (;;)
    {
    char *e;
    if (s == NULL || s[0] == 0 || count == maxArraySize)
	break;
    e = strchr(s, ',');
    if (e != NULL)
	*e++ = 0;
    array[count++] = atof(s);
    s = e;
    }
return count;
}

void sqlFloatDynamicArray(char *s, float **retArray, int *retSize)
/* Convert comma separated list of numbers to an dynamically allocated
 * array, which should be freeMem()'d when done. Thread-safe. */
{
float *array = NULL;
int count = 0;

if (s)
    {
    count = countSeparatedItems(s, ',');
    if (count > 0)
	{
	AllocArray(array, count);
	count = 0;
	for (;;)
	    {
	    array[count++] = sqlFloatInList(&s);
	    if (*s++ == 0)
		break;
	    if (*s == 0)
		break;
	    }
	}
    }
*retArray = array;
*retSize = count;
}

void sqlSignedDynamicArray(char *s, int **retArray, int *retSize)
/* Convert comma separated list of numbers to an dynamically allocated
 * array, which should be freeMem()'d when done. Thread-safe. */
{
int *array = NULL;
int count = 0;

if (s)
    {
    count = countSeparatedItems(s, ',');
    if (count > 0)
	{
	AllocArray(array, count);
	count = 0;
	for (;;)
	    {
	    array[count++] = sqlSignedInList(&s);
	    if (*s++ == 0)
		break;
	    if (*s == 0)
		break;
	    }
	}
    }
*retArray = array;
*retSize = count;
}

int sqlUnsignedComma(char **pS)
/* Return signed number at *pS.  Advance *pS past comma at end */
{
char *s = *pS;
char *e = strchr(s, ',');
unsigned ret;

*e++ = 0;
*pS = e;
ret = sqlUnsigned(s);
return ret;
}


int sqlSignedComma(char **pS)
/* Return signed number at *pS.  Advance *pS past comma at end */
{
char *s = *pS;
char *e = strchr(s, ',');
int ret;

*e++ = 0;
*pS = e;
ret = sqlSigned(s);
return ret;
}

float sqlFloatComma(char **pS)
/* Return signed number at *pS.  Advance *pS past comma at end */
{
char *s = *pS;
char *e = strchr(s, ',');
float ret;

*e++ = 0;
*pS = e;
ret = atof(s);
return ret;
}

static char *findStringEnd(char *start, char endC)
/* Return end of string. */
{
char c;
char *s = start;

for (;;)
    {
    c = *s;
    if (c == endC)
	return s;
    else if (c == 0)
	errAbort("Unterminated string");
    ++s;
    }
}

static char *sqlGetOptQuoteString(char **pS)
/* Return string at *pS.  (Either quoted or not.)  Advance *pS. */
{
char *s = *pS;
char *e;
char c = *s;

if (c  == '"' || c == '\'')
    {
    s += 1;
    e = findStringEnd(s, c);
    *e++ = 0;
    if (*e++ != ',')
	errAbort("Expecting comma after string");
    }
else
    {
    e = strchr(s, ',');
    *e++ = 0;
    }
*pS = e;
return s;
}

char *sqlStringComma(char **pS)
/* Return string at *pS.  (Either quoted or not.)  Advance *pS. */
{
return cloneString(sqlGetOptQuoteString(pS));
}

void sqlFixedStringComma(char **pS, char *buf, int bufSize)
/* Copy string at *pS to buf.  Advance *pS. */
{
strncpy(buf, sqlGetOptQuoteString(pS), bufSize);
}

char *sqlEatChar(char *s, char c)
/* Make sure next character is 'c'.  Return past next char */
{
if (*s++ != c)
    errAbort("Expecting %c got %c (%d) in database", c, s[-1], s[-1]);
return s;
}

static struct hash *buildSymHash(char **values, boolean isEnum)
/* build a hash of values for either enum or set symbolic column */
{
struct hash *valHash = hashNew(0);
unsigned setVal = 1; /* not used for enum */
int iVal;
for (iVal = 0; values[iVal] != NULL; iVal++)
    {
    if (isEnum)
        hashAddInt(valHash, values[iVal], iVal);
    else
        {
        hashAddInt(valHash, values[iVal], setVal);
        setVal = setVal << 1;
        }
    }
return valHash;
}

unsigned sqlEnumParse(char *valStr, char **values, struct hash **valHashPtr)
/* parse an enumerated column value */
{
if (*valHashPtr == NULL)
    *valHashPtr = buildSymHash(values, TRUE);
return hashIntVal(*valHashPtr, valStr);
}

unsigned sqlSetParse(char *valStr, char **values, struct hash **valHashPtr)
/* parse a set column value */
{
if (*valHashPtr == NULL)
    *valHashPtr = buildSymHash(values, FALSE);
/* parse comma separated string */
unsigned value = 0;
char *val = strtok(valStr, ",");
while (val != NULL)
    {
    value |= hashIntVal(*valHashPtr, val);
    val = strtok(NULL, ",");
    }

return value;
}
