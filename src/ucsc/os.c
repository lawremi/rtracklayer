/* Added by rtracklayer to condition on the platform */

/* some shared utilities */

#include "common.h"
#include "portable.h"

struct fileInfo *newFileInfo(char *name, off_t size, bool isDir, int statErrno, 
                             time_t lastAccess)
/* Return a new fileInfo. */
{
  int len = strlen(name);
  struct fileInfo *fi = needMem(sizeof(*fi) + len);
  fi->size = size;
  fi->isDir = isDir;
  fi->statErrno = statErrno;
  fi->lastAccess = lastAccess;
  strcpy(fi->name, name);
  return fi;
}

int cmpFileInfo(const void *va, const void *vb)
/* Compare two fileInfo. */
{
  const struct fileInfo *a = *((struct fileInfo **)va);
  const struct fileInfo *b = *((struct fileInfo **)vb);
  return strcmp(a->name, b->name);
}


#ifdef WIN32
#include "oswin9x.c"
#else
#include "osunix.c"
#endif
