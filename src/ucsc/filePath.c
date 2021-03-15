 /* filePath - stuff to handle file name parsing. */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "filePath.h"

static char *findSlashBefore(char *start, char *e)
/* Return first slash before s (but not before start) */
{
    while (--e >= start)
    {
        if (*e == '/')
            return e;
    }
    return start;
}

void undosPath(char *path)
/* Convert '\' to '/' in path. */
{
    subChar(path, '\\', '/');
}

char *expandRelativePath(char *baseDir, char *relPath)
/* Expand relative path to more absolute one. */
{
    if (relPath[0] == '/')
    // hey, it's absolute actually... 
        return cloneString(relPath);

    char *e = baseDir + strlen(baseDir);
    int slashCount;
    char *rel = relPath;
    char *result;
    int size, baseSize;
    undosPath(baseDir);
    undosPath(relPath);
    slashCount = countChars(baseDir, '/');
    if (baseDir[0] == 0)
        slashCount = -1;
    while (startsWith("../", rel))
        {
            if (slashCount < 0)
            {
                warn("More ..'s in \"%s\" than directories in \"%s\"", relPath, baseDir);
                return NULL;
            }
            else if (slashCount == 0)
                e = baseDir;
            else
                e = findSlashBefore(baseDir, e);
        slashCount -= 1;
        rel += 3;
        }
    baseSize = e - baseDir;
    size = strlen(rel) + 1;
    if (baseSize > 0)
        size += baseSize + 1;
    if (baseSize > 0)
        {
            result = needMem(size);
            memcpy(result, baseDir, baseSize);
            result[baseSize] = '/';
            strcpy(result + baseSize + 1, rel);
        }
    else
        result = cloneString(rel);
    return result;
}
