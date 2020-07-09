/* memgfx - routines for drawing on bitmaps in memory.
 * Currently limited to 256 color bitmaps. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "memgfx.h"


struct rgbColor colorIxToRgb(int colorIx)
/* Return rgb value at color index. */
{
static struct rgbColor rgb;
#ifdef MEMGFX_BIGENDIAN
rgb.r = (colorIx >> 24) & 0xff;
rgb.g = (colorIx >> 16) & 0xff;
rgb.b = (colorIx >> 8) & 0xff;
#else
rgb.r = (colorIx >> 0) & 0xff;
rgb.g = (colorIx >> 8) & 0xff;
rgb.b = (colorIx >> 16) & 0xff;
#endif
return rgb;
}

