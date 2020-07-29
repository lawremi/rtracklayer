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

  SEXP ans, n_qhits, ranges, chromStart, chromWidth, name, score, strand,
    thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts,
    expCount, expIds, expScores, label, extraFields, lengthIndex, extraNames;

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
  int fieldCount = file->fieldCount;
  int definedFieldCount = file->definedFieldCount;
  int extraFieldCount = fieldCount - definedFieldCount;
  int n_hits = slCount(hits);
  SEXPTYPE *typeId;
  char *asText = bigBedAutoSqlText(file);
  struct asObject *as = asParseText(asText);
  freeMem(asText);
  bigBedFileClose(&file);

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
      }
      asCol = asCol->next;
    }
    lengthIndex = PROTECT(allocVector(INTSXP, extraFieldCount));
    memset(INTEGER(lengthIndex), 0, sizeof(int) * extraFieldCount);
  }

  chromStart = PROTECT(allocVector(INTSXP, n_hits));
  chromWidth = PROTECT(allocVector(INTSXP, n_hits));
  name = PROTECT(allocVector(STRSXP, n_hits));
  score = PROTECT(allocVector(INTSXP, n_hits));
  strand = PROTECT(allocVector(STRSXP, n_hits));
  thickStart = PROTECT(allocVector(INTSXP, n_hits));
  thickEnd = PROTECT(allocVector(INTSXP, n_hits));
  itemRgb = PROTECT(allocVector(STRSXP, n_hits));
  blockCount = PROTECT(allocVector(INTSXP, n_hits));
  blockSizes = PROTECT(allocVector(INTSXP, n_hits));
  blockStarts = PROTECT(allocVector(INTSXP, n_hits));

  char startBuf[16], endBuf[16], *row[fieldCount], rgbBuf[8];
  for (int i = 0, k = 0; i < n_hits; ++i, hits = hits->next) {
    if (INTEGER(n_qhits)[k] == i && k < n_ranges)
      ++k;
    bigBedIntervalToRow(hits, (char *)CHAR(STRING_ELT(r_seqnames, k)),
                        startBuf, endBuf, row, fieldCount);
    struct bed *bed = bedLoadN(row, definedFieldCount);
    /* storing default fields */
    INTEGER(chromStart)[i] = bed->chromStart;
    INTEGER(chromWidth)[i] = bed->chromEnd - bed->chromStart + 1;
    SET_STRING_ELT(name, i, mkChar(bed->name));
    INTEGER(score)[i] = bed->score;
    SET_STRING_ELT(strand, i, mkChar(bed->strand));
    INTEGER(thickStart)[i] = bed->thickStart;
    INTEGER(thickEnd)[i] = bed->thickEnd;
    snprintf(rgbBuf, 8, "#%x", bed->itemRgb);
    SET_STRING_ELT(itemRgb, i, mkChar(rgbBuf));
    INTEGER(blockCount)[i] = bed->blockCount;
    INTEGER(blockSizes)[i] = bed->blockSizes;
    INTEGER(blockStarts)[i] = bed->chromStarts;
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
  ans = PROTECT(allocVector(VECSXP, 13));

  SET_VECTOR_ELT(ans, 0, ranges);
  SET_VECTOR_ELT(ans, 1, name);
  SET_VECTOR_ELT(ans, 2, score);
  SET_VECTOR_ELT(ans, 3, strand);
  SET_VECTOR_ELT(ans, 4, thickStart);
  SET_VECTOR_ELT(ans, 5, thickEnd);
  SET_VECTOR_ELT(ans, 6, itemRgb);
  SET_VECTOR_ELT(ans, 7, blockCount);
  SET_VECTOR_ELT(ans, 8, blockSizes);
  SET_VECTOR_ELT(ans, 9, blockStarts);
  SET_VECTOR_ELT(ans, 10, extraFields);
  SET_VECTOR_ELT(ans, 11, extraNames);
  SET_VECTOR_ELT(ans, 12, n_qhits);

  extraFieldCount > 0 ? UNPROTECT(17 + extraFieldCount) : UNPROTECT(16);
  lmCleanup(&lm);
  popRHandlers();
  return ans;
}
