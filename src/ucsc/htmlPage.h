#ifndef HTMLPAGE_H
#define HTMLPAGE_H

#ifndef DYSTRING_H
#include "dystring.h"
#endif

char *expandUrlOnBase(char *base, char *url);
/* Figure out first character past host name. Load up
 * return string with protocol (if any) and host name. 
 * It is assumed that url is relative to base and does not contain a protocol.*/

#endif /* HTMLPAGE_H */
