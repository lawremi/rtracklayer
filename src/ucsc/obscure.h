/* Obscure.h  - stuff that's relatively rarely used
 * but still handy. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef OBSCURE_H
#define OBSCURE_H

int digitsBaseTwo(unsigned long x);
/* Return base two # of digits. */

void sprintLongWithCommas(char *s, long long l);
/* Print out a long number with commas a thousands, millions, etc. */

void writeGulp(char *file, char *buf, int size);
/* Write out a bunch of memory. */

void readInGulp(char *fileName, char **retBuf, size_t *retSize);
/* Read whole file in one big gulp. */

void *intToPt(int i);
/* Convert integer to pointer. Use when really want to store an
 * int in a pointer field. */

int ptToInt(void *pt);
/* Convert pointer to integer.  Use when really want to store a
 * pointer in an int. */

boolean parseQuotedString( char *in, char *out, char **retNext);
/* Read quoted string from in (which should begin with first quote).
 * Write unquoted string to out, which may be the same as in.
 * Return pointer to character past end of string in *retNext. 
 * Return FALSE if can't find end. */

void escCopy(char *in, char *out, char toEscape, char escape);
/* Copy in to out, escaping as needed.  Out better be big enough. 
 * (Worst case is strlen(in)*2 + 1.) */

struct slName *charSepToSlNames(char *string, char c);
/* Convert character-separated list of items to slName list. 
 * Note that the last occurence of c is optional.  (That
 * is for a comma-separated list a,b,c and a,b,c, are
 * equivalent. */

struct hash *hashThisEqThatLine(char *line, int lineIx, boolean firstStartsWithLetter);
/* Return a symbol table from a line of form:
 *   1-this1=val1 2-this='quoted val2' var3="another val" 
 * If firstStartsWithLetter is true, then the left side of the equals must start with
 * and equals. */


void shuffleArrayOfPointers(void *pointerArray, int arraySize);
/* Shuffle array of pointers of given size given number of times. */

void shuffleList(void *pList);
/* Randomize order of slList.  Usage:
 *     shuffleList(&list)
 * where list is a pointer to a structure that
 * begins with a next field. */

void *slListRandomReduce(void *list, double reduceRatio);
/* Reduce list to approximately reduceRatio times original size. Destroys original list. */

void rangeRoundUp(double start, double end, double *retStart, double *retEnd);
/* Round start and end so that they cover a slightly bigger range, but with more round
 * numbers.  For instance 0.23:9.89 becomes 0:10 */

void rangeFromMinMaxMeanStd(double minVal, double maxVal, double mean, double std,
	double *retStart, double *retEnd);
/* Given some basic statistical properties, set a range that will be good on a wide
 * range of biological data. */

#endif /* OBSCURE_H */
