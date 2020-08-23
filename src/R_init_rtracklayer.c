#include "rtracklayer.h"
#include "readGFF.h"
#include "bigWig.h"
#include "bigBed.h"
#include "twoBit.h"
#include "utils.h"

#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* readGFF.c */
  CALLMETHOD_DEF(gff_colnames, 1),
  CALLMETHOD_DEF(read_gff_pragmas, 1),
  CALLMETHOD_DEF(scan_gff, 5),
  CALLMETHOD_DEF(load_gff, 8),
  /* bigWig.c */
  CALLMETHOD_DEF(BWGSectionList_add, 5),
  CALLMETHOD_DEF(BWGSectionList_write, 5),
  CALLMETHOD_DEF(BWGSectionList_cleanup, 1),
  CALLMETHOD_DEF(BWGFile_query, 5),
  CALLMETHOD_DEF(BWGFile_seqlengths, 1),
  CALLMETHOD_DEF(BWGFile_summary, 6),
  CALLMETHOD_DEF(BWGFile_fromWIG, 4),
  CALLMETHOD_DEF(R_udcCleanup, 1),
  CALLMETHOD_DEF(R_setUserUdcDir, 1),
  /* bigBed.c */
  CALLMETHOD_DEF(BBDFile_fieldnames, 1),
  CALLMETHOD_DEF(BBDFile_seqlengths, 1),
  CALLMETHOD_DEF(BBDFile_query, 5),
  CALLMETHOD_DEF(BBDFile_write, 6),
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
