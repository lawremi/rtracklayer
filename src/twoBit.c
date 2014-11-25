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

/* .Call entry point */
/* Writes the list of twoBit pointers to disk and frees them */
SEXP TwoBits_write(SEXP r_twoBits, SEXP r_filename) {
  pushRHandlers();
  FILE *file = mustOpen((char *)CHAR(asChar(r_filename)), "wb");
  struct twoBit *twoBits = NULL, *twoBit_it = NULL;
  
  for (int i = 0; i < length(r_twoBits); i++)
    slAddHead(&twoBits, R_ExternalPtrAddr(VECTOR_ELT(r_twoBits, i)));
  slReverse(&twoBits);
  
  twoBitWriteHeader(twoBits, file);
  for (twoBit_it = twoBits; twoBit_it != NULL; twoBit_it = twoBit_it->next) {
    twoBitWriteOne(twoBit_it, file);
  }
  
  twoBitFreeList(&twoBits);
  carefulClose(&file);
  popRHandlers();
  
  return R_NilValue;
}

/* .Call entry point */
SEXP TwoBitFile_seqlengths(SEXP r_filename) {
  pushRHandlers();
  struct twoBitFile *tbf = twoBitOpen((char *)CHAR(asChar(r_filename)));
  struct twoBitIndex *index;
  int i, n = slCount(tbf->indexList);
  SEXP r_seqlengths, r_seqnames;
  
  PROTECT(r_seqlengths = allocVector(INTSXP, n));
  r_seqnames = allocVector(STRSXP, n);
  setAttrib(r_seqlengths, R_NamesSymbol, r_seqnames);
  
  for (index = tbf->indexList, i = 0; index != NULL; index = index->next, i++) {
    SET_STRING_ELT(r_seqnames, i, mkChar(index->name));
    INTEGER(r_seqlengths)[i] = twoBitSeqSize(tbf, index->name);
  }

  twoBitClose(&tbf);
  popRHandlers();
  UNPROTECT(1);
  
  return r_seqlengths;
}

SEXP TwoBitFile_read(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges, SEXP lkup)
{
  pushRHandlers();
  struct twoBitFile *file = twoBitOpen((char *)CHAR(asChar(r_filename)));
  int *frag_start = INTEGER(get_IRanges_start(r_ranges));
  int *frag_width = INTEGER(get_IRanges_width(r_ranges));
  int frag_count = get_IRanges_length(r_ranges);
  SEXP r_ans_width, r_ans;
  XVectorList_holder r_ans_holder;

  PROTECT(r_ans_width = duplicate(get_IRanges_width(r_ranges)));
  PROTECT(r_ans = alloc_XRawList("DNAStringSet", "DNAString", r_ans_width));
  r_ans_holder = hold_XVectorList(r_ans);
  for (int i = 0; i < frag_count; i++) {
    if (frag_width[i]) { // UCSC library does not like zero width ranges
      struct dnaSeq *frag =
        twoBitReadSeqFrag(file, (char *)CHAR(STRING_ELT(r_seqnames, i)),
                          frag_start[i] - 1, frag_start[i] + frag_width[i] - 1);
      Chars_holder r_ans_elt_holder =
        get_elt_from_XRawList_holder(&r_ans_holder, i);
      /* r_ans_elt_holder.ptr is a const char * so we need to cast it to
         char * before we can write to it */
      Ocopy_bytes_to_i1i2_with_lkup(0, r_ans_elt_holder.length - 1,
        (char *)r_ans_elt_holder.ptr, r_ans_elt_holder.length,
        frag->dna, frag->size,
        INTEGER(lkup), LENGTH(lkup));
      freeDnaSeq(&frag);
    }
  }
  twoBitClose(&file);
  popRHandlers();

  UNPROTECT(2);
  return r_ans;
}
