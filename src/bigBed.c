#include "ucsc/common.h"
#include "ucsc/linefile.h"
#include "ucsc/localmem.h"
#include "ucsc/bbiFile.h"
#include "ucsc/bigBed.h"
#include "ucsc/basicBed.h"
#include "ucsc/asParse.h"
#include "ucsc/sqlNum.h"

#include "bigBed.h"
#include "bbiHelper.h"
#include "handlers.h"

enum IFields
{
  i_name    = 3,  /* Index value of name field    */
  i_score   = 4,  /* Index value of score field   */
  i_strand  = 5,  /* Index value of strand field  */
  i_thick   = 7,  /* Index value of thick field   */
  i_itemRgb = 8,  /* Index value of itemRgb field */
  i_blocks  = 11, /* Index value of blocks field  */
};

/* helper functions */
int getDefinedFieldCount(struct asObject *as);
bool isPresent(int definedFieldCount, int index);

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_seqlengths(SEXP r_filename)
{
  pushRHandlers();
  struct bbiFile *file = bigBedFileOpen((char *)CHAR(asChar(r_filename)));
  SEXP seqlengths = PROTECT(bbiSeqLengths(file));
  bigBedFileClose(&file);
  popRHandlers();
  UNPROTECT(1);
  return seqlengths;
}

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_query(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges)
{
  pushRHandlers();
  struct bbiFile *file = bigBedFileOpen((char *)CHAR(asChar(r_filename)));
  struct lm *lm = lmInit(0);
  int n_ranges = get_IRanges_length(r_ranges);
  int *start = INTEGER(get_IRanges_start(r_ranges));
  int *width = INTEGER(get_IRanges_width(r_ranges));

  SEXP ans, n_qhits, ranges, chromStart, chromWidth, name, score,
    strand = R_NilValue, thickStart, thickEnd, itemRgb, blockCount, blockSizes,
    blockStarts, expCount, expIds, expScores, label, extraFields, lengthIndex,
    extraNames;

  n_qhits = PROTECT(allocVector(INTSXP, n_ranges));
  struct bigBedInterval *hits = NULL, *tail = NULL;
  /* querying records in range */
  for (int i = 0; i < n_ranges; ++i) {
    struct bigBedInterval *queryHits =
      bigBedIntervalQuery(file, (char *)CHAR(STRING_ELT(r_seqnames, i)),
                          start[i] - 1, start[i] - 1 + width[i], 0, lm);
    if (!hits) {
      hits = queryHits;
      tail = slLastEl(hits);
    } else {
      tail->next = queryHits;
      tail = slLastEl(tail);
    }
    INTEGER(n_qhits)[i] = slCount(queryHits);
  }

  /* need these before closing file */
  char *asText = bigBedAutoSqlText(file);
  struct asObject *as = asParseText(asText);
  freeMem(asText);
  int fieldCount = file->fieldCount;
  int definedFieldCount = getDefinedFieldCount(as);
  int extraFieldCount = fieldCount - definedFieldCount;
  int n_hits = slCount(hits);
  bigBedFileClose(&file);

  int presentFieldCount = 0, unprotectCount = 0;
  /* mandatory default field */
  chromStart = PROTECT(allocVector(INTSXP, n_hits));
  chromWidth = PROTECT(allocVector(INTSXP, n_hits));
  /* if any of the default field is present allocate memory for it */
  if (isPresent(definedFieldCount, i_name)) {
    name = PROTECT(allocVector(STRSXP, n_hits));
    ++presentFieldCount;
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_score)) {
    score = PROTECT(allocVector(INTSXP, n_hits));
    ++presentFieldCount;
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_strand)) {
    strand = PROTECT(allocVector(STRSXP, n_hits));
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_thick)) {
    thickStart = PROTECT(allocVector(INTSXP, n_hits));
    thickEnd = PROTECT(allocVector(INTSXP, n_hits));
    presentFieldCount += 2;
    unprotectCount += 2;
  }
  if (isPresent(definedFieldCount, i_itemRgb)) {
    itemRgb = PROTECT(allocVector(STRSXP, n_hits));
    ++presentFieldCount;
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_blocks)) {
      blockCount = PROTECT(allocVector(INTSXP, n_hits));
      blockSizes = PROTECT(allocVector(INTSXP, n_hits));
      blockStarts = PROTECT(allocVector(INTSXP, n_hits));
      presentFieldCount += 3;
      unprotectCount += 3;
  }

  SEXPTYPE *typeId;
  extraFields = PROTECT(allocVector(VECSXP, extraFieldCount));
  extraNames = PROTECT(allocVector(STRSXP, extraFieldCount));

  /* storing extra field type and extra names */
  if (extraFieldCount > 0) {
    typeId = (SEXPTYPE*)R_alloc(extraFieldCount, sizeof(SEXPTYPE));
    struct asColumn *asCol = as->columnList;
    enum asTypes fieldType;
    for (int j = 0, k = 0; j < fieldCount; ++j) {
      fieldType = asCol->lowType->type;
      if (j >= definedFieldCount) {
        SET_STRING_ELT(extraNames, k, mkChar(asCol->name));
        if (fieldType >= 0 && fieldType <= 1)
          typeId[k] = REALSXP;
        else if (fieldType >= 3 && fieldType <= 9)
          typeId[k] = INTSXP;
        else if ((fieldType == 2) || (fieldType >= 10 && fieldType <= 15))
          typeId[k] = STRSXP;
        SEXP temp = PROTECT(allocVector(typeId[k], n_hits));
        SET_VECTOR_ELT(extraFields, k, temp);
        ++k;
        ++unprotectCount;
      }
      asCol = asCol->next;
    }
    lengthIndex = PROTECT(allocVector(INTSXP, extraFieldCount));
    memset(INTEGER(lengthIndex), 0, sizeof(int) * extraFieldCount);
    ++unprotectCount;
  }
  asObjectFree(&as);

  char startBuf[16], endBuf[16], *row[fieldCount], rgbBuf[8];
  for (int i = 0, k = 0; i < n_hits; ++i, hits = hits->next) {
    if (INTEGER(n_qhits)[k] == i && k < n_ranges)
      ++k;
    bigBedIntervalToRow(hits, (char *)CHAR(STRING_ELT(r_seqnames, k)),
                        startBuf, endBuf, row, fieldCount);
    struct bed *bed = bedLoadN(row, definedFieldCount);
    /* mandatory default field */
    INTEGER(chromStart)[i] = bed->chromStart;
    INTEGER(chromWidth)[i] = bed->chromEnd - bed->chromStart + 1;
    /* if any of the default field is present, store its values */
    if (isPresent(definedFieldCount, i_name)) {
      SET_STRING_ELT(name, i, mkChar(bed->name));
    }
    if (isPresent(definedFieldCount, i_score)) {
      INTEGER(score)[i] = bed->score;
    }
    if (isPresent(definedFieldCount, i_strand)) {
      SET_STRING_ELT(strand, i, mkChar(bed->strand));
    }
    if (isPresent(definedFieldCount, i_thick)) {
      INTEGER(thickStart)[i] = bed->thickStart;
      INTEGER(thickEnd)[i] = bed->thickEnd;
    }
    if (isPresent(definedFieldCount, i_itemRgb)) {
      snprintf(rgbBuf, 8, "#%x", bed->itemRgb);
      SET_STRING_ELT(itemRgb, i, mkChar(rgbBuf));
    }
    if (isPresent(definedFieldCount, i_blocks)) {
      INTEGER(blockCount)[i] = bed->blockCount;
      INTEGER(blockSizes)[i] = bed->blockSizes;
      INTEGER(blockStarts)[i] = bed->chromStarts;
    }
    bedFree(&bed);

    /* if extra fields are present store extra fields */
    for (int j = definedFieldCount, efIndex = 0 ; j < fieldCount; ++j, ++efIndex) {
      switch(typeId[efIndex]) {
        case REALSXP:
          REAL(VECTOR_ELT(extraFields, efIndex))[i] = sqlSigned(row[j]);
          break;
        case INTSXP:
          INTEGER(VECTOR_ELT(extraFields, efIndex))[i] = sqlUnsigned(row[j]);
          break;
        case STRSXP: {
          int index = INTEGER(lengthIndex)[efIndex];
          SET_STRING_ELT(VECTOR_ELT(extraFields, efIndex), index, mkChar(row[j]));
          INTEGER(lengthIndex)[efIndex] = index + 1;
          break;
        }
      }
    }
  }

  ranges = PROTECT(new_IRanges("IRanges", chromStart, chromWidth, R_NilValue));
  ans = PROTECT(allocVector(VECSXP, presentFieldCount + 5));
  int index = 0;
  SET_VECTOR_ELT(ans, index++, n_qhits);
  SET_VECTOR_ELT(ans, index++, extraFields);
  SET_VECTOR_ELT(ans, index++, extraNames);
  SET_VECTOR_ELT(ans, index++, ranges);
  SET_VECTOR_ELT(ans, index++, strand);
  if (isPresent(definedFieldCount, i_name)) {
    SET_VECTOR_ELT(ans, index++, name);
  }
  if (isPresent(definedFieldCount, i_score)) {
    SET_VECTOR_ELT(ans, index++, score);
  }
  if (isPresent(definedFieldCount, i_thick)) {
    SET_VECTOR_ELT(ans, index++, thickStart);
    SET_VECTOR_ELT(ans, index++, thickEnd);
  }
  if (isPresent(definedFieldCount, i_itemRgb)) {
    SET_VECTOR_ELT(ans, index++, itemRgb);
  }
  if (isPresent(definedFieldCount, i_blocks)) {
    SET_VECTOR_ELT(ans, index++, blockCount);
    SET_VECTOR_ELT(ans, index++, blockSizes);
    SET_VECTOR_ELT(ans, index++, blockStarts);
  }
  UNPROTECT(7 + unprotectCount);
  lmCleanup(&lm);
  popRHandlers();
  return ans;
}

int getDefinedFieldCount(struct asObject *as) {
  int definedFieldCount = 0;
  struct asColumn *asCol = as->columnList;
  char *asText = bedAsDef(12, 12);
  struct asObject *bedAs = asParseText(asText);
  freeMem(asText);
  struct asColumn *bedCol = bedAs->columnList;
  while (asCol) {
    if (strncmp(asCol->name, bedCol->name, sizeof(asCol->name)) == 0)
      ++definedFieldCount;
    bedCol = bedCol->next;
    asCol = asCol->next;
  }
  asObjectFree(&bedAs);
  return definedFieldCount;
}

bool isPresent(int definedFieldCount, int index) {
  return (definedFieldCount - index >= 0) ? TRUE : FALSE;
}
