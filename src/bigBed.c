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
bool isSelected(SEXP r_selectedindex, int position);

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
SEXP BBDFile_fieldnames(SEXP r_filename)
{
  pushRHandlers();
  struct bbiFile *file = bigBedFileOpen((char *)CHAR(asChar(r_filename)));
  char *asText = bigBedAutoSqlText(file);
  struct asObject *as = asParseText(asText);
  freeMem(asText);
  int fieldCount = file->fieldCount;
  int definedFieldCount = getDefinedFieldCount(as);
  bigBedFileClose(&file);
  char *names[] = {"name", "score", "thick", "itemRgb", "blocks"};
  struct asColumn *asCol = as->columnList;
  SEXP defaultFields = PROTECT(allocVector(STRSXP, definedFieldCount));
  SEXP extraFields = PROTECT(allocVector(STRSXP, fieldCount - definedFieldCount));
  for (int i = 0; i < fieldCount; ++i) {
    if (i >= definedFieldCount)
      SET_STRING_ELT(extraFields, i - definedFieldCount, mkChar(asCol->name));
    else if (i == 3)
      SET_STRING_ELT(defaultFields, i, mkChar(names[0]));
    else if (i == 4)
      SET_STRING_ELT(defaultFields, i, mkChar(names[1]));
    else if (i == 7)
      SET_STRING_ELT(defaultFields, i, mkChar(names[2]));
    else if (i == 8)
      SET_STRING_ELT(defaultFields, i, mkChar(names[3]));
    else if (i == 11)
      SET_STRING_ELT(defaultFields, i, mkChar(names[4]));
    asCol = asCol->next;
  }
  SEXP list = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(list, 0, defaultFields);
  SET_VECTOR_ELT(list, 1, extraFields);
  asObjectFree(&as);
  popRHandlers();
  UNPROTECT(3);
  return list;
}

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_query(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges,
                   SEXP r_defaultindex, SEXP r_extraindex)
{
  pushRHandlers();
  struct bbiFile *file = bigBedFileOpen((char *)CHAR(asChar(r_filename)));
  struct lm *lm = lmInit(0);
  int n_ranges = get_IRanges_length(r_ranges);
  int *start = INTEGER(get_IRanges_start(r_ranges));
  int *width = INTEGER(get_IRanges_width(r_ranges));

  SEXP ans, n_qhits, ranges, chromStart, chromWidth, name, score,
    strand = R_NilValue, thickStart, thickWidth, itemRgb, blocks,
    extraFields = R_NilValue, lengthIndex;

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
  /* if any of the default field is present and selected, allocate memory for it */
  if (isPresent(definedFieldCount, i_name) && isSelected(r_defaultindex, 1)) {
    name = PROTECT(allocVector(STRSXP, n_hits));
    ++presentFieldCount;
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_score) && isSelected(r_defaultindex, 2)) {
    score = PROTECT(allocVector(INTSXP, n_hits));
    ++presentFieldCount;
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_strand)) {
    strand = PROTECT(allocVector(STRSXP, n_hits));
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_thick) && isSelected(r_defaultindex, 3)) {
    thickStart = PROTECT(allocVector(INTSXP, n_hits));
    thickWidth = PROTECT(allocVector(INTSXP, n_hits));
    ++presentFieldCount;
    unprotectCount += 2;
  }
  if (isPresent(definedFieldCount, i_itemRgb) && isSelected(r_defaultindex, 4)) {
    itemRgb = PROTECT(allocVector(STRSXP, n_hits));
    ++presentFieldCount;
    ++unprotectCount;
  }
  if (isPresent(definedFieldCount, i_blocks) && isSelected(r_defaultindex, 5)) {
      blocks = PROTECT(allocVector(VECSXP, n_hits));
      ++presentFieldCount;
      ++unprotectCount;
  }

  SEXPTYPE *typeId;
  /* if extra fields are present and selected
   * identify the type information and allocate memory */
  if (extraFieldCount > 0) {
    int k = 0, i = 0;
    enum asTypes fieldType;
    struct asColumn *asCol = as->columnList;
    extraFields = PROTECT(allocVector(VECSXP, extraFieldCount));
    typeId = (SEXPTYPE*)R_alloc(extraFieldCount, sizeof(SEXPTYPE));
    for (int j = 0; j < fieldCount; ++j) {
      fieldType = asCol->lowType->type;
      if (j >= definedFieldCount) {
        if (asTypesIsFloating(fieldType) || fieldType == t_int ||
            fieldType == t_short || fieldType == t_byte || fieldType == t_off) {
          typeId[k] = REALSXP;
        } else if (asTypesIsUnsigned(fieldType)) {
          typeId[k] = INTSXP;
        } else if (fieldType == t_char || fieldType == t_string ||
                   fieldType == t_lstring) {
          typeId[k] = STRSXP;
        }
        if (isSelected(r_extraindex, (j - definedFieldCount + 1))) {
          SEXP temp = PROTECT(allocVector(typeId[k], n_hits));
          SET_VECTOR_ELT(extraFields, k, temp);
          ++unprotectCount;
          ++k;
        }
      }
      asCol = asCol->next;
    }
    lengthIndex = PROTECT(allocVector(INTSXP, extraFieldCount));
    memset(INTEGER(lengthIndex), 0, sizeof(int) * extraFieldCount);
    unprotectCount += 2;
  }
  asObjectFree(&as);

  int count = 0, k = 0;
  char startBuf[16], endBuf[16], *row[fieldCount], rgbBuf[8];
  for (int i = 0; i < n_hits; ++i, hits = hits->next, ++count) {
    if (INTEGER(n_qhits)[k] == count && k < n_ranges) {
      ++k;
      count = 0;
    }
    bigBedIntervalToRow(hits, (char *)CHAR(STRING_ELT(r_seqnames, k)),
                        startBuf, endBuf, row, fieldCount);
    struct bed *bed = bedLoadN(row, definedFieldCount);
    /* mandatory default field */
    INTEGER(chromStart)[i] = bed->chromStart;
    INTEGER(chromWidth)[i] = bed->chromEnd - bed->chromStart + 1;
    /* if any of the default field is present and selected, store its value */
    if (isPresent(definedFieldCount, i_name) && isSelected(r_defaultindex, 1)) {
      SET_STRING_ELT(name, i, mkChar(bed->name));
    }
    if (isPresent(definedFieldCount, i_score) && isSelected(r_defaultindex, 2)) {
      INTEGER(score)[i] = bed->score;
    }
    if (isPresent(definedFieldCount, i_strand)) {
      SET_STRING_ELT(strand, i, mkChar(bed->strand));
    }
    if (isPresent(definedFieldCount, i_thick) && isSelected(r_defaultindex, 3)) {
      INTEGER(thickWidth)[i] = bed->thickEnd - bed->thickStart + 1;
      INTEGER(thickStart)[i] = bed->thickStart;
    }
    if (isPresent(definedFieldCount, i_itemRgb) && isSelected(r_defaultindex, 4)) {
      snprintf(rgbBuf, 8, "#%x", bed->itemRgb);
      SET_STRING_ELT(itemRgb, i, mkChar(rgbBuf));
    }
    if (isPresent(definedFieldCount, i_blocks) && isSelected(r_defaultindex, 5)) {
      SEXP bstart = PROTECT(allocVector(INTSXP, bed->blockCount));
      SEXP bwidth = PROTECT(allocVector(INTSXP, bed->blockCount));
      for (int j = 0; j< bed->blockCount; ++j) {
        INTEGER(bwidth)[j] = bed->blockSizes[j];
        INTEGER(bstart)[j] = bed->chromStarts[j];
      }
      SET_VECTOR_ELT(blocks, i, new_IRanges("IRanges", bstart, bwidth, R_NilValue));
      UNPROTECT(2);
    }
    bedFree(&bed);

    /* if extra fields are present and selected store their values */
    for (int j = definedFieldCount, efIndex = 0 ; j < fieldCount; ++j) {
      if (isSelected(r_extraindex, (j - definedFieldCount + 1))) {
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
        ++efIndex;
      }
    }
    freeMem(row[3]);
  }

  ranges = PROTECT(new_IRanges("IRanges", chromStart, chromWidth, R_NilValue));
  ans = PROTECT(allocVector(VECSXP, presentFieldCount + 4));
  int index = 0;
  SET_VECTOR_ELT(ans, index++, n_qhits);
  SET_VECTOR_ELT(ans, index++, extraFields);
  SET_VECTOR_ELT(ans, index++, ranges);
  SET_VECTOR_ELT(ans, index++, strand);
  if (isPresent(definedFieldCount, i_name) && isSelected(r_defaultindex, 1)) {
    SET_VECTOR_ELT(ans, index++, name);
  }
  if (isPresent(definedFieldCount, i_score) && isSelected(r_defaultindex, 2)) {
    SET_VECTOR_ELT(ans, index++, score);
  }
  if (isPresent(definedFieldCount, i_thick) && isSelected(r_defaultindex, 3)) {
    SET_VECTOR_ELT(ans, index++, new_IRanges("IRanges", thickStart,
                                             thickWidth, R_NilValue));
  }
  if (isPresent(definedFieldCount, i_itemRgb) && isSelected(r_defaultindex, 4)) {
    SET_VECTOR_ELT(ans, index++, itemRgb);
  }
  if (isPresent(definedFieldCount, i_blocks) && isSelected(r_defaultindex, 5)) {
    SET_VECTOR_ELT(ans, index++, blocks);
  }
  UNPROTECT(5 + unprotectCount);
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
  while (asCol && bedCol) {
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

bool isSelected(SEXP r_selectedindex, int position) {
  if (length(r_selectedindex) == 0)
    return TRUE;
  for (int i = 0; i < length(r_selectedindex); ++i) {
    if(INTEGER(r_selectedindex)[i] == position)
      return TRUE;
  }
  return FALSE;
}
