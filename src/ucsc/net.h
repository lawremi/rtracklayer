/* Net.h some stuff to wrap around net communications. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */


#ifndef NET_H
#define NET_H

#include "linefile.h"
#include "dystring.h"

#define DEFAULTCONNECTTIMEOUTMSEC 10000  /* default connect timeout for tcp in milliseconds */
#define DEFAULTREADWRITETTIMEOUTSEC 120  /* default read/write timeout for tcp in seconds */
#define MAXURLSIZE 4096 /* maximum size in characters for a URL, but also see the struct netParsedUrl definition */

int netUrlHead(char *url, struct hash *hash, char **effectiveUrl);
/* Go get head and return status. Returns HTTP status code, or 0 on other errors.
 * If hash is non-null, fill it with header lines. Keys are lower-cased.
 * If effectiveUrl is non-NULL, it will be set to a clone of the final URL
 * after redirects if a redirect occurred, otherwise it will be NULL.
 * The caller is responsible for freeing effectiveUrl. */

int netUrlRead(char *url, bits64 offset, int size, void *buffer);
/* Fetch a range of bytes from a URL.
 * Returns bytes read, or -1 on error. */

boolean netGetFtpInfo(char *url, long long *retSize, time_t *retTime);
/* Return date in UTC and size of ftp url file */

boolean hasProtocol(char *urlOrPath);
/* Return TRUE if it looks like it has http://, ftp:// etc. */
#endif /* NET_H */

