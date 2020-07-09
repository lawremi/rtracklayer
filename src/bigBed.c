#include "ucsc/common.h"
#include "ucsc/linefile.h"
#include "ucsc/localmem.h"
#include "ucsc/bbiFile.h"
#include "ucsc/bigBed.h"

#include "bigBed.h"
#include "handlers.h"

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_seqlengths(SEXP r_filename) {
  pushRHandlers();
  struct bbiFile * file = bigBedFileOpen((char *)CHAR(asChar(r_filename)));
  struct bbiChromInfo *chromList = bbiChromList(file);
  struct bbiChromInfo *chrom = chromList;
  SEXP seqlengths, seqlengthNames;

  PROTECT(seqlengths = allocVector(INTSXP, slCount(chromList)));
  seqlengthNames = allocVector(STRSXP, length(seqlengths));
  setAttrib(seqlengths, R_NamesSymbol, seqlengthNames);

  for(int i = 0; i < length(seqlengths); i++) {
    INTEGER(seqlengths)[i] = chrom->size;
    SET_STRING_ELT(seqlengthNames, i, mkChar(chrom->name));
    chrom = chrom->next;
  }

  bigBedFileClose(&file);
  bbiChromInfoFreeList(&chromList);
  popRHandlers();
  UNPROTECT(1);
  return seqlengths;
}
