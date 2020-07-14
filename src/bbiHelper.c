#include "ucsc/common.h"
#include "ucsc/bbiFile.h"

#include "bbiHelper.h"

SEXP bbiSeqLengths(struct bbiFile *file) {
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
  bbiChromInfoFreeList(&chromList);
  UNPROTECT(1);
  return seqlengths;
}
