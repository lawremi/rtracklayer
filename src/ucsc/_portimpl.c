/* rtracklayer took this utility out of portimpl.c */

#include "common.h"
#include "portable.h"

void makeDirsOnPath(char *pathName)
/* Create directory specified by pathName.  If pathName contains
 * slashes, create directory at each level of path if it doesn't
 * already exist.  Abort with error message if there's a problem.
 * (It's not considered a problem for the directory to already
 * exist. ) */
{
/* Save a copy of current directory. */
char *curDir = cloneString(getCurrentDir());

/* Return current directory.  Abort if it fails. */
/* Make local copy of pathName. */
int len = strlen(pathName);
char pathCopy[len+1];
strcpy(pathCopy, pathName);

/* Start at root if it's an absolute path name. */
char *s = pathCopy, *e;
if (pathCopy[0] == '/')
  {
setCurrentDir("/");
++s;
}

/* Step through it one slash at a time - changing directory if possible
 * else making directory if possible, else dying. */
for (; !isEmpty(s); s = e)
  {
/* Find end of this section and terminate string there. */
e = strchr(s, '/');
if (e != NULL)
  *e++ = 0;

/* Cd there.  If that fails mkdir there and then cd there. */
if (!maybeSetCurrentDir(s))
  {
if (!makeDir(s))
  {
break;
}
setCurrentDir(s);
}
}
setCurrentDir(curDir);
freeMem(curDir);
}
