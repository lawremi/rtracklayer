#include "rtracklayer.h"
#include "bigWig.h"
#include "twoBit.h"
#include "utils.h"

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
  CALLMETHOD_DEF(BWGFile_fromWIG, 3),
  /* twobit.c */
  CALLMETHOD_DEF(DNAString_to_twoBit, 3),
  CALLMETHOD_DEF(TwoBits_write, 2),
  CALLMETHOD_DEF(TwoBitFile_seqlengths, 1),
  CALLMETHOD_DEF(TwoBitFile_read, 4),
  /* utils.c */
  CALLMETHOD_DEF(CharacterList_pasteCollapse, 2),
  {NULL, NULL, 0}
};

void R_init_rtracklayer(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
