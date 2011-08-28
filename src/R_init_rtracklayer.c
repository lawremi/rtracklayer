#include "rtracklayer.h"
#include "bigWig.h"

#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
/* bigWig.c */
  CALLMETHOD_DEF(BWGSectionList_add, 5),
  CALLMETHOD_DEF(BWGSectionList_write, 4),
  CALLMETHOD_DEF(BWGSectionList_cleanup, 1),
  CALLMETHOD_DEF(BWGFile_query, 3),
  CALLMETHOD_DEF(BWGFile_seqlengths, 1),
  CALLMETHOD_DEF(BWGFile_summary, 6),
  {NULL, NULL, 0}
};

void R_init_rtracklayer(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
