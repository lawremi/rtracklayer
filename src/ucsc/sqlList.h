/* Stuff for processing comma separated lists .
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef SQLLIST_H
#define SQLLIST_H
struct hash;

int sqlDoubleArray(char *s, double *array, int maxArraySize);
/* Return sum of double values in a comma-separated list */

int sqlFloatArray(char *s, float *array, int maxArraySize);

void sqlFloatDynamicArray(char *s, float **retArray, int *retSize);
void sqlSignedDynamicArray(char *s, int **retArray, int *retSize);

char *sqlEscapeString(const char *orig);
/* Prepares string for inclusion in a SQL statement . Remember to free
 * returned string.  returned string contains strlen(length)*2+1 as many bytes
 * as orig because in worst case every character has to be escaped.
 * Example 1: The Gene's Name -> The Gene''s Name
 * Example 2: he said "order and orient" -> he said ""order and orient"" */

char *sqlEscapeString2(char *to, const char* from);
/* Prepares a string for inclusion in a sql statement.  Output string
 * must be 2*strlen(from)+1 */

int sqlUnsignedComma(char **pS);
/* Return signed number at *pS.  Advance *pS past comma at end.
 * This function is used by the system that automatically loads
 * structured object from longblobs. */

int sqlSignedComma(char **pS);
/* Return signed number at *pS.  Advance *pS past comma at end */

float sqlFloatComma(char **pS);
/* Return floating point number at *pS.  Advance *pS past comma at end */

char *sqlStringComma(char **pS);
/* Return string at *pS.  (Either quoted or not.)  Advance *pS. */

void sqlFixedStringComma(char **pS, char *buf, int bufSize);
/* Copy string at *pS to buf.  Advance *pS. */

char *sqlEatChar(char *s, char c);
/* Make sure next character is 'c'.  Return past next char */

unsigned sqlEnumParse(char *valStr, char **values, struct hash **valHashPtr);
/* parse an enumerated column value */

unsigned sqlSetParse(char *valStr, char **values, struct hash **valHashPtr);
/* parse a set column value */

#endif /* SQLLIST_H */

