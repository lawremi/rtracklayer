/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "hash.h"
#include "dystring.h"
#include "filePath.h"
#include "htmlPage.h"

char *expandUrlOnBase(char *base, char *url)
/* Figure out first character past host name. Load up
 * return string with protocol (if any) and host name. 
 * It is assumed that url is relative to base and does not contain a protocol.*/
{
    struct dyString *dy = NULL;
    char *hostName, *pastHostName;
    dy = dyStringNew(256);
    if (startsWith("http:", base) || startsWith("https:", base) || startsWith("ftp:", base))
        hostName = (strchr(base, ':') + 3);
    else
        hostName = base;
    pastHostName = strchr(hostName, '/');
    if (pastHostName == NULL)
        pastHostName = hostName + strlen(hostName);
    dyStringAppendN(dy, base, pastHostName - base);

    /* Add url to return string after host name. */
    if (startsWith("/", url))	/* New URL is absolute, just append to hostName */
    {
        dyStringAppend(dy, url);
    }
    else
    {
        char *curDir = pastHostName;
        char *endDir;
        if (curDir[0] == '/')
            curDir += 1;
        dyStringAppendC(dy, '/');
        endDir = strrchr(curDir, '/');
        if (endDir == NULL)
            endDir = curDir;
        if (startsWith("../", url))
        {
            char *dir = cloneStringZ(curDir, endDir-curDir);
            char *path = expandRelativePath(dir, url);
            if (path != NULL)
            {
                dyStringAppend(dy, path);
            }
            freez(&dir);
            freez(&path);
        }
        else
        {
            dyStringAppendN(dy, curDir, endDir-curDir);
            if (lastChar(dy->string) != '/')
                dyStringAppendC(dy, '/');
            dyStringAppend(dy, url);
        }
    }
    return dyStringCannibalize(&dy);
}
