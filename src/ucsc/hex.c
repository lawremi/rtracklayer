/* Handy hexidecimal functions
 *   If you don't want to use printf
 */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"

char hexTab[16] = {'0', '1', '2', '3', '4', '5', '6', '7', 
	'8', '9', 'a', 'b', 'c', 'd', 'e', 'f', };
/* Convert 0-15 to a hex char */

void hexBinaryString(unsigned char *in, int inSize, char *out, int outSize)
/* Convert possibly long binary string to hex string.
 * Out size needs to be at least 2x inSize+1 */
{
    assert(inSize * 2 +1 <= outSize);
    while (--inSize >= 0)
    {
        unsigned char c = *in++;
        *out++ = hexTab[c>>4];
        *out++ = hexTab[c&0xf];
    }
    *out = 0;
}
