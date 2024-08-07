#include "ucsc/common.h"
#include "ucsc/dnaseq.h"
#include "ucsc/twoBit.h"

#include "twoBit.h"
#include "handlers.h"

/* .Call entry point */
SEXP DNAString_to_twoBit(SEXP r_dna, SEXP r_mask, SEXP r_seqname) {
  pushRHandlers();
  dnaUtilOpen();
  const DNA *dna = CHAR(asChar(r_dna));
  struct dnaSeq *seq = newDnaSeq((DNA *)dna, strlen(dna),
                                 (char *)CHAR(asChar(r_seqname)));
  struct twoBit *twoBit = twoBitFromDnaSeq(seq, FALSE);
  int *mask_start = INTEGER(get_IRanges_start(r_mask));
  int *mask_width = INTEGER(get_IRanges_width(r_mask));
  int mask_count = get_IRanges_length(r_mask);
  SEXP ans;

  if (mask_count) {
    AllocArray(twoBit->maskStarts, mask_count);
    AllocArray(twoBit->maskSizes, mask_count);
  }
  for (int i = 0; i < mask_count; i++) {
    twoBit->maskStarts[i] = mask_start[i] - 1;
    twoBit->maskSizes[i] = mask_width[i];
  }

  seq->dna = NULL; /* do not free memory owned by R */
  freeDnaSeq(&seq);
  popRHandlers();

  PROTECT(ans = R_MakeExternalPtr(twoBit, R_NilValue, R_NilValue));
  setAttrib(ans, R_ClassSymbol, mkString("twoBit"));
  UNPROTECT(1);

  return ans;
}
