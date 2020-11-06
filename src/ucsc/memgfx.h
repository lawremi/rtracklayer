/* Memgfx - stuff to do graphics in memory buffers.
 * Typically will just write these out as .gif or .png files.
 * This stuff is byte-a-pixel for simplicity.
 * It can do 256 colors.
 *
 * This file is copyright 2000 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef MEMGFX_H
#define MEMGFX_H

#if defined(__sgi__) || defined(__sgi) || defined(__powerpc__) || defined(sparc) || defined(__ppc__) || defined(__s390__) || defined(__s390x__)

// BIGENDIAN machines:

#define MEMGFX_BIGENDIAN	1
#define MG_WHITE   0xffffffff
#define MG_BLACK   0x000000ff
#define MG_RED     0xff0000ff
#define MG_GREEN   0x00ff00ff
#define MG_BLUE    0x0000ffff
#define MG_CYAN    0x00ffffff
#define MG_MAGENTA 0xff00ffff
#define MG_YELLOW  0xffff00ff
#define MG_GRAY    0x808080ff

#define MAKECOLOR_32(r,g,b) (((unsigned int)0xff) | ((unsigned int)b<<8) | ((unsigned int)g << 16) | ((unsigned int)r << 24))
#define COLOR_32_RED(c) (((c)>>24)&0xff)
#define COLOR_32_GREEN(c) (((c)>>16)&0xff)
#define COLOR_32_BLUE(c) (((c)>>8)&0xff)

#else

// LITTLE ENDIAN machines:

#define MG_WHITE   0xffffffff
#define MG_BLACK   0xff000000
#define MG_RED     0xff0000ff
#define MG_GREEN   0xff00ff00
#define MG_BLUE    0xffff0000
#define MG_CYAN    0xffffff00
#define MG_MAGENTA 0xffff00ff
#define MG_YELLOW  0xff00ffff
#define MG_GRAY    0xff808080

#define MAKECOLOR_32(r,g,b) (((unsigned int)0xff<<24) | ((unsigned int)b<<16) | ((unsigned int)g << 8) | (unsigned int)r)
#define COLOR_32_RED(c) ((c)&0xff)
#define COLOR_32_GREEN(c) (((c)>>8)&0xff)
#define COLOR_32_BLUE(c) (((c)>>16)&0xff)
#endif

struct rgbColor
    {
    unsigned char r, g, b;
    };

struct rgbColor colorIxToRgb(int colorIx);
/* Return rgb value at color index. */

#endif /* MEMGFX_H */
