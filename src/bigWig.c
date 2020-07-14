#include "ucsc/common.h"
#include "ucsc/linefile.h"
#include "ucsc/localmem.h"
#include "ucsc/hash.h"
#include "ucsc/bbiFile.h"
#include "ucsc/bigWig.h"
#include "ucsc/bwgInternal.h"

#include "bigWig.h"
#include "bbiHelper.h"
#include "handlers.h"

static struct bwgBedGraphItem *
createBedGraphItems(int *start, int *width, double *score, int len,
                    struct lm *lm)
{
  struct bwgBedGraphItem *itemList = NULL, *item;
  int i;
  for (i=0; i<len; ++i)
    {
      lmAllocVar(lm, item);
      item->end = start[i] + width[i] - 1;
      item->start = start[i] - 1;
      item->val = score[i];
      slAddHead(&itemList, item);
    }
  slReverse(&itemList);
  return itemList;
}

static struct bwgVariableStepPacked *
createVariableStepItems(int *start, double *score, int len, struct lm *lm)
{
  struct bwgVariableStepPacked *packed;
  lmAllocArray(lm, packed, len);
  int i;
  for (i=0; i<len; ++i)
    {
      packed[i].start = start[i] - 1;
      packed[i].val = score[i];
    }
  return packed;
}

static struct bwgFixedStepPacked *
createFixedStepItems(double *score, int len, struct lm *lm)
{
  struct bwgFixedStepPacked *packed;
  lmAllocArray(lm, packed, len);
  int i;
  for (i=0; i<len; ++i)
    {
      packed[i].val = score[i];
    }
  return packed;
}

static struct bwgSection *
createBWGSection_Rle(const char *seq, int *start, int *width, double *score,
                     int len, enum bwgSectionType type, struct lm *lm)
{
  struct bwgSection *section;
  lmAllocVar(lm, section);
  section->chrom = (char *)seq;
  section->start = start[0] - 1;
  section->end = start[len-1] + width[len-1] - 1;
  section->type = type;
  section->itemSpan = width[0];
  if (type == bwgTypeFixedStep) {
    section->items.fixedStepPacked = createFixedStepItems(score, len, lm);
    section->itemStep = len > 1 ? start[1] - start[0] : 0;
  } else if (type == bwgTypeVariableStep) {
    section->items.variableStepPacked =
      createVariableStepItems(start, score, len, lm);
  } else section->items.bedGraphList =
           createBedGraphItems(start, width, score, len, lm);
  section->itemCount = len;
  return section;
}

static struct bwgSection *
createBWGSection_Atomic(const char *seq, int start, double *score,
                        int len, struct lm *lm)
{
  struct bwgSection *section;
  lmAllocVar(lm, section);
  section->chrom = (char *)seq;
  section->start = start;
  section->end = start + len - 1;
  section->type = bwgTypeFixedStep;
  section->itemSpan = 1;
  section->items.fixedStepPacked = createFixedStepItems(score, len, lm);
  section->itemStep = 1;
  section->itemCount = len;
  return section;
}

static int itemsPerSlot = 1024;
static int blockSize = 256;

static void BWGSectionList_addRle(struct bwgSection **sections, const char *seq,
                                  SEXP r_ranges, double *score,
                                  enum bwgSectionType type, struct lm *lm)
{
  int numLeft = get_IRanges_length(r_ranges);
  int *start = INTEGER(get_IRanges_start(r_ranges));
  int *width = INTEGER(get_IRanges_width(r_ranges));
  while(numLeft) {
    int numSection = numLeft > itemsPerSlot ? itemsPerSlot : numLeft;
    numLeft -= numSection;
    slAddHead(sections, createBWGSection_Rle(seq, start, width, score,
                                             numSection, type, lm));
    start += numSection;
    width += numSection;
    score += numSection;
  }
}

static void BWGSectionList_addAtomic(struct bwgSection **sections,
                                     const char *seq,
                                     double *score,
                                     int num,
                                     struct lm *lm)
{
  int numLeft = num;
  while(numLeft) {
    int numSection = numLeft > itemsPerSlot ? itemsPerSlot : numLeft;
    slAddHead(sections, createBWGSection_Atomic(seq, num - numLeft, score,
                                                numSection, lm));
    score += numSection;
    numLeft -= numSection;
  }
}

/* --- .Call ENTRY POINT --- */

SEXP BWGSectionList_add(SEXP r_sections, SEXP r_seq, SEXP r_ranges,
                        SEXP r_score, SEXP r_format)
{
  struct bwgSection *sections = NULL;
  const char *seq = CHAR(asChar(r_seq));
  double *score = REAL(r_score);
  const char *format = CHAR(asChar(r_format));
  SEXP ans;
  struct lm *lm;

  enum bwgSectionType type = bwgTypeBedGraph;
  if (sameString(format, "fixedStep"))
    type = bwgTypeFixedStep;
  else if (sameString(format, "variableStep"))
    type = bwgTypeVariableStep;

  if (r_sections != R_NilValue) {
    sections = R_ExternalPtrAddr(r_sections);
    lm = R_ExternalPtrAddr(R_ExternalPtrTag(r_sections));
  } else lm = lmInit(0);

  pushRHandlers();
  if (r_ranges != R_NilValue) {
    BWGSectionList_addRle(&sections, seq, r_ranges, score, type, lm);
  } else {
    BWGSectionList_addAtomic(&sections, seq, score, length(r_score), lm);
  }
  popRHandlers();

  PROTECT(ans = R_MakeExternalPtr(sections, R_NilValue, R_NilValue));
  R_SetExternalPtrTag(ans, R_MakeExternalPtr(lm, R_NilValue, R_NilValue));
  UNPROTECT(1);

  return ans;
}

static struct hash *createIntHash(SEXP v) {
  struct hash *hash = hashNew(0);
  SEXP names = getAttrib(v, R_NamesSymbol);
  for (int i = 0; i < length(v); i++)
    hashAddInt(hash, (char *)CHAR(STRING_ELT(names, i)), INTEGER(v)[i]); 
  return hash;
}

/* --- .Call ENTRY POINT --- */
SEXP BWGSectionList_write(SEXP r_sections, SEXP r_seqlengths, SEXP r_compress,
                          SEXP r_fixed_summaries, SEXP r_file)
{
  struct bwgSection *sections = NULL;
  struct hash *lenHash = createIntHash(r_seqlengths);
  if (r_sections != R_NilValue) {
    sections = R_ExternalPtrAddr(r_sections);
    slReverse(&sections);
  }
  pushRHandlers();
  bwgCreate(sections, lenHash, max(blockSize, length(r_seqlengths)),
	    itemsPerSlot, asLogical(r_compress),
            FALSE /*keepAllChromosomes*/, asLogical(r_fixed_summaries),
            (char *)CHAR(asChar(r_file)));
  freeHash(&lenHash);
  popRHandlers();
  return r_file;
}

/* --- .Call ENTRY POINT --- */
SEXP BWGSectionList_cleanup(SEXP r_sections)
{
  pushRHandlers();
  if (r_sections != R_NilValue) {
    struct lm *lm = R_ExternalPtrAddr(R_ExternalPtrTag(r_sections));
    lmCleanup(&lm);
  }
  popRHandlers();
  return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP BWGFile_seqlengths(SEXP r_filename) {
  pushRHandlers();
  SEXP seqlengths;
  struct bbiFile * file = bigWigFileOpen((char *)CHAR(asChar(r_filename)));
  PROTECT(seqlengths = bbiSeqLengths(file));
  bbiFileClose(&file); 
  popRHandlers();
  UNPROTECT(1);
  return seqlengths;
}

/* --- .Call ENTRY POINT --- */
SEXP BWGFile_query(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges,
		   SEXP r_return_score, SEXP r_return_list) {
  pushRHandlers();
  struct bbiFile * file = bigWigFileOpen((char *)CHAR(asChar(r_filename)));
  Rboolean return_list = asLogical(r_return_list);
  SEXP ans, ans_start, ans_width, ans_score, ans_ranges, ans_nhits;
  SEXP numericListEls = NULL;
  bool returnScore = asLogical(r_return_score);
  struct lm *lm = lmInit(0);
  
  struct bbiInterval *hits = NULL;

  int n_ranges = get_IRanges_length(r_ranges);

  if (return_list) {
      PROTECT(numericListEls = allocVector(VECSXP, n_ranges));
  } else {
      PROTECT(ans_nhits = allocVector(INTSXP, n_ranges));
  }
  
  int *start = INTEGER(get_IRanges_start(r_ranges));
  int *width = INTEGER(get_IRanges_width(r_ranges));

  for (int i = 0; i < n_ranges; i++) {
      struct bbiInterval *queryHits =
	  bigWigIntervalQuery(file, (char *)CHAR(STRING_ELT(r_seqnames, i)),
			      start[i] - 1, start[i] - 1 + width[i], lm);
      int nqhits = slCount(queryHits);
      if (return_list) {
	  struct bbiInterval *qhits = queryHits;
	  SEXP ans_numeric;
	  PROTECT(ans_numeric = allocVector(REALSXP, width[i]));
	  memset(REAL(ans_numeric), 0, sizeof(double) * width[i]);
	  for (int k = 0; k < nqhits; k++, qhits = qhits->next) {
	      for (int l = qhits->start; l < qhits->end; l++)
		  REAL(ans_numeric)[(l - start[i] + 1)] = qhits->val;
	  }
	  SET_VECTOR_ELT(numericListEls, i, ans_numeric);
	  UNPROTECT(1);
      } else {
	  slReverse(&queryHits);
	  hits = slCat(queryHits, hits);
	  INTEGER(ans_nhits)[i] = nqhits;
      }
  }

  bbiFileClose(&file);
  
  if (return_list) {
    ans = new_SimpleList("SimpleList", numericListEls);
    UNPROTECT(1);
  } else {
    int n_hits = slCount(hits);
    PROTECT(ans_start = allocVector(INTSXP, n_hits));
    PROTECT(ans_width = allocVector(INTSXP, n_hits));
    if (returnScore) {
      PROTECT(ans_score = allocVector(REALSXP, n_hits));
    } else ans_score = R_NilValue;

    slReverse(&hits);
    for (int i = 0; i < n_hits; i++, hits = hits->next) {
      INTEGER(ans_start)[i] = hits->start + 1;
      INTEGER(ans_width)[i] = hits->end - hits->start;
      if (returnScore)
	REAL(ans_score)[i] = hits->val;
    }

    PROTECT(ans_ranges = new_IRanges("IRanges", ans_start, ans_width,
	R_NilValue));

    ans = allocVector(VECSXP, 3);
    SET_ELEMENT(ans, 0, ans_ranges);
    SET_ELEMENT(ans, 1, ans_score);
    SET_ELEMENT(ans, 2, ans_nhits);

    UNPROTECT(4 + returnScore);
  }

  lmCleanup(&lm);
  popRHandlers();
  return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP BWGFile_summary(SEXP r_filename, SEXP r_chrom, SEXP r_ranges,
                     SEXP r_size, SEXP r_type, SEXP r_default_value)
{
  pushRHandlers();
  struct bbiFile * file = bigWigFileOpen((char *)CHAR(asChar(r_filename)));
  enum bbiSummaryType type =
    bbiSummaryTypeFromString((char *)CHAR(asChar(r_type)));
  double default_value = asReal(r_default_value);
  int *start = INTEGER(get_IRanges_start(r_ranges));
  int *width = INTEGER(get_IRanges_width(r_ranges));
  SEXP ans;
  
  PROTECT(ans = allocVector(VECSXP, length(r_chrom)));
  for (int i = 0; i < length(r_chrom); i++) {
    int size = INTEGER(r_size)[i];
    char *chrom = (char *)CHAR(STRING_ELT(r_chrom, i));
    SEXP r_values = allocVector(REALSXP, size);
    double *values = REAL(r_values);
    for (int j = 0; j < size; j++)
      values[j] = default_value;
    SET_VECTOR_ELT(ans, i, r_values);
    bool success = bigWigSummaryArray(file, chrom, start[i] - 1,
                                      start[i] - 1 + width[i], type, size,
                                      values);
    if (!success)
      warning("Failed to summarize range %d (%s:%d-%d)", i, chrom, start[i],
            start[i] - 1 + width[i]);
  }
  bbiFileClose(&file);
  popRHandlers();
  UNPROTECT(1);
  return ans;
}

#include "ucsc/verbose.h"

SEXP BWGFile_fromWIG(SEXP r_infile, SEXP r_clip, SEXP r_seqlengths,
		     SEXP r_outfile)
{
  pushRHandlers();
  struct lm *lm = lmInit(0);
  struct hash *lenHash = createIntHash(r_seqlengths);
  Rboolean clip = asLogical(r_clip);
  struct bwgSection *sections =
    bwgParseWig((char *)CHAR(asChar(r_infile)), clip, lenHash, itemsPerSlot,
                lm);
  bwgCreate(sections, lenHash, max(blockSize, length(r_seqlengths)),
	    itemsPerSlot, TRUE, TRUE, FALSE,
            (char *)CHAR(asChar(r_outfile)));
  lmCleanup(&lm);
  freeHash(&lenHash);
  popRHandlers();
  return r_outfile;
}

SEXP R_udcCleanup(SEXP r_maxDays) {
    double maxDays = asReal(r_maxDays);
    bits64 size = udcCleanup(udcDefaultDir(), maxDays, FALSE);
    return ScalarReal(size);
}

SEXP R_setUserUdcDir(SEXP r_dir) {
  char *dir = (char *)CHAR(asChar(r_dir));
  udcSetDefaultDir(dir);
  return R_NilValue;
}
