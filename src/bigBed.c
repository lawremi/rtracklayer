#include "ucsc/common.h"
#include "ucsc/hash.h"
#include "ucsc/bigBed.h"
#include "ucsc/linefile.h"
#include "ucsc/localmem.h"

#include "bigBed.h"
#include "handlers.h"
#include "bbiHelper.h"
#include "bigBedHelper.h"

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
    int k = 0;
    enum asTypes fieldType;
    struct asColumn *asCol = as->columnList;
    extraFields = PROTECT(allocVector(VECSXP, extraFieldCount));
    typeId = (SEXPTYPE*)R_alloc(extraFieldCount, sizeof(SEXPTYPE));
    for (int j = 0; j < fieldCount; ++j) {
      fieldType = asCol->lowType->type;
      if (j >= definedFieldCount) {
        if (asTypesIsFloating(fieldType) || fieldType == t_uint ||
            fieldType == t_off) {
          typeId[k] = REALSXP;
        } else if (fieldType == t_int || fieldType == t_short ||
                   fieldType == t_ushort || fieldType == t_byte) {
          typeId[k] = INTSXP;
        } else if (fieldType == t_char || fieldType == t_string ||
                   fieldType == t_lstring) {
          typeId[k] = STRSXP;
        } else if (fieldType == t_ubyte) {
          typeId[k] = RAWSXP;
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
      snprintf(rgbBuf, 8, "#%06x", bed->itemRgb);
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
            REAL(VECTOR_ELT(extraFields, efIndex))[i] = sqlDouble(row[j]);
            break;
          case INTSXP:
            INTEGER(VECTOR_ELT(extraFields, efIndex))[i] = sqlSigned(row[j]);
            break;
          case STRSXP: {
            int index = INTEGER(lengthIndex)[efIndex];
            SET_STRING_ELT(VECTOR_ELT(extraFields, efIndex), index, mkChar(row[j]));
            INTEGER(lengthIndex)[efIndex] = index + 1;
            break;
          }
          case RAWSXP:
            RAW(extraFields)[i] = sqlUnsigned(row[j]);
            break;
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

static struct hash *createIntHash(SEXP v) {
  struct hash *hash = hashNew(0);
  SEXP names = getAttrib(v, R_NamesSymbol);
  for (int i = 0; i < length(v); ++i)
    hashAddInt(hash, (char *)CHAR(STRING_ELT(names, i)), INTEGER(v)[i]);
  return hash;
}

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_write(SEXP r_seqlengths, SEXP r_bedString, SEXP r_autosql,
                   SEXP r_indexfields, SEXP r_compress, SEXP r_outfile)
{
  pushRHandlers();
  int blockSize = 256;
  int itemsPerSlot = 512;
  char *bedString = cloneString((char *)CHAR(asChar(r_bedString)));
  struct lineFile *lf = lineFileOnString("text", TRUE, bedString);
  struct bbExIndexMaker *eim = NULL;
  bool doCompress = asLogical(r_compress);
  struct hash *lenHash = createIntHash(r_seqlengths);
  char *asText = (char *)CHAR(asChar(r_autosql));
  struct asObject *as = asParseText(asText);
  bits16 fieldCount = slCount(as->columnList);
  bits16 definedFieldCount = getDefinedFieldCount(as);
  char *extraIndex = (char *)CHAR(asChar(r_indexfields));
  struct slName *extraIndexList = slNameListFromString(extraIndex, ',');
  bits16 extraIndexCount = slCount(extraIndexList);
  if (extraIndexList != NULL)
    eim = bbExIndexMakerNew(extraIndexList, as);

  /* Do first pass, mostly just scanning file and counting hits per chromosome. */
  int minDiff = 0;
  double aveSize = 0;
  bits64 bedCount = 0;
  bits32 uncompressBufSize = 0;
  struct bbiChromUsage *usageList = bbiChromUsageFromBedFile(lf, lenHash, eim, &minDiff,
                                                             &aveSize, &bedCount);

  /* Open output file and write dummy header. */
  FILE *f = mustOpen((char *)CHAR(asChar(r_outfile)), "wb");
  bbiWriteDummyHeader(f);
  bbiWriteDummyZooms(f);

  /* Write out autoSql string */
  bits64 asOffset = ftell(f);
  mustWrite(f, asText, strlen(asText) + 1);

  /* Write out dummy total summary. */
  struct bbiSummaryElement totalSum;
  ZeroVar(&totalSum);
  bits64 totalSummaryOffset = ftell(f);
  bbiSummaryElementWrite(f, &totalSum);

  /* Write out dummy header extension */
  bits64 extHeaderOffset = ftell(f);
  bits16 extHeaderSize = 64;
  repeatCharOut(f, 0, extHeaderSize);

  /* Write out extra index stuff if need be. */
  bits64 extraIndexListOffset = 0;
  bits64 extraIndexListEndOffset = 0;
  if (extraIndexList != NULL) {
    extraIndexListOffset = ftell(f);
    int extraIndexSize = 16 + 4*1;   /* Fixed record size 16, plus 1 times field size of 4 */
    repeatCharOut(f, 0, extraIndexSize*extraIndexCount);
    extraIndexListEndOffset = ftell(f);
  }

  /* Write out chromosome/size database. */
  bits64 chromTreeOffset = ftell(f);
  bbiWriteChromInfo(usageList, blockSize, f);

  /* Set up to keep track of possible initial reduction levels. */
  int resScales[bbiMaxZoomLevels], resSizes[bbiMaxZoomLevels];
  int resTryCount = bbiCalcResScalesAndSizes(aveSize, resScales, resSizes);

  /* Write out primary full resolution data in sections, collect stats to use for reductions. */
  bits64 dataOffset = ftell(f);
  bits32 blockCount = 0;
  bits32 maxBlockSize = 0;
  struct bbiBoundsArray *boundsArray = NULL;
  writeOne(f, bedCount);
  if (bedCount > 0) {
    blockCount = bbiCountSectionsNeeded(usageList, itemsPerSlot);
    AllocArray(boundsArray, blockCount);
    freez(&bedString);
    bedString = cloneString((char *)CHAR(asChar(r_bedString)));
    lf = lineFileOnString("text", TRUE, bedString);
    if (eim)
      bbExIndexMakerAllocChunkArrays(eim, bedCount);
    writeBlocks(usageList, lf, as, itemsPerSlot, boundsArray, blockCount, doCompress,
                f, resTryCount, resScales, resSizes, eim, bedCount, fieldCount,
                definedFieldCount, &maxBlockSize);
  }

  /* Write out primary data index. */
  bits64 indexOffset = ftell(f);
  cirTreeFileBulkIndexToOpenFile(boundsArray, sizeof(boundsArray[0]), blockCount,
                                 blockSize, 1, NULL, bbiBoundsArrayFetchKey,
                                 bbiBoundsArrayFetchOffset, indexOffset, f);
  freez(&boundsArray);

  /* Declare arrays and vars that track the zoom levels we actually output. */
  bits32 zoomAmounts[bbiMaxZoomLevels];
  bits64 zoomDataOffsets[bbiMaxZoomLevels];
  bits64 zoomIndexOffsets[bbiMaxZoomLevels];

  /* Call monster zoom maker library function that bedGraphToBigWig also uses. */
  int zoomLevels = 0;
  if (bedCount > 0) {
    freez(&bedString);
    bedString = cloneString((char *)CHAR(asChar(r_bedString)));
    lf = lineFileOnString("text", TRUE, bedString);
    zoomLevels = bbiWriteZoomLevels(lf, f, blockSize, itemsPerSlot, bedWriteReducedOnceReturnReducedTwice,
                                    fieldCount, doCompress, indexOffset - dataOffset, usageList,
                                    resTryCount, resScales, resSizes, zoomAmounts, zoomDataOffsets,
                                    zoomIndexOffsets, &totalSum);
  }

  /* Write out extra indexes if need be. */
  if (eim) {
    int i;
    for (i=0; i < eim->indexCount; ++i) {
      eim->fileOffsets[i] = ftell(f);
      maxBedNameSize = eim->maxFieldSize[i];
      qsort(eim->chunkArrayArray[i], bedCount,
            sizeof(struct bbNamedFileChunk), bbNamedFileChunkCmpByName);
      assert(sizeof(struct bbNamedFileChunk) == sizeof(eim->chunkArrayArray[i][0]));
      bptFileBulkIndexToOpenFile(eim->chunkArrayArray[i], sizeof(eim->chunkArrayArray[i][0]),
                                 bedCount, blockSize, bbNamedFileChunkKey, maxBedNameSize,
                                 bbNamedFileChunkVal, sizeof(bits64) + sizeof(bits64), f);
    }
  }

  /* Figure out buffer size needed for uncompression if need be. */
  if (doCompress) {
    int maxZoomUncompSize = itemsPerSlot * sizeof(struct bbiSummaryOnDisk);
    uncompressBufSize = max(maxBlockSize, maxZoomUncompSize);
  }

  /* Go back and rewrite header. */
  rewind(f);
  bits32 sig = bigBedSig;
  bits16 version = bbiCurrentVersion;
  bits16 summaryCount = zoomLevels;
  bits32 reserved32 = 0;
  bits64 reserved64 = 0;

  /* Write fixed header */
  writeOne(f, sig);
  writeOne(f, version);
  writeOne(f, summaryCount);
  writeOne(f, chromTreeOffset);
  writeOne(f, dataOffset);
  writeOne(f, indexOffset);
  writeOne(f, fieldCount);
  writeOne(f, definedFieldCount);
  writeOne(f, asOffset);
  writeOne(f, totalSummaryOffset);
  writeOne(f, uncompressBufSize);
  writeOne(f, extHeaderOffset);
  assert(ftell(f) == 64);

  /* Write summary headers with data. */
  int i;
  for (i=0; i<zoomLevels; ++i) {
    writeOne(f, zoomAmounts[i]);
    writeOne(f, reserved32);
    writeOne(f, zoomDataOffsets[i]);
    writeOne(f, zoomIndexOffsets[i]);
  }
  /* Write rest of summary headers with no data. */
  for (i=zoomLevels; i<bbiMaxZoomLevels; ++i) {
    writeOne(f, reserved32);
    writeOne(f, reserved32);
    writeOne(f, reserved64);
    writeOne(f, reserved64);
  }

  /* Write total summary. */
  fseek(f, totalSummaryOffset, SEEK_SET);
  bbiSummaryElementWrite(f, &totalSum);

  /* Write extended header */
  fseek(f, extHeaderOffset, SEEK_SET);
  writeOne(f, extHeaderSize);
  writeOne(f, extraIndexCount);
  writeOne(f, extraIndexListOffset);
  repeatCharOut(f, 0, 52);    // reserved
  assert(ftell(f) - extHeaderOffset == extHeaderSize);

  /* Write extra index offsets if need be. */
  if (extraIndexCount != 0) {
    fseek(f, extraIndexListOffset, SEEK_SET);
    int i;
    for (i = 0; i < extraIndexCount; ++i) {
      // Write out fixed part of index info
      bits16 type = 0;    // bPlusTree type
      bits16 indexFieldCount = 1;
      writeOne(f, type);
      writeOne(f, indexFieldCount);
      writeOne(f, eim->fileOffsets[i]);
      repeatCharOut(f, 0, 4);  // reserved

      // Write out field list - easy this time because for now always only one field.
      bits16 fieldId = eim->indexFields[i];
      writeOne(f, fieldId);
      repeatCharOut(f, 0, 2); // reserved
    }
    assert(ftell(f) == extraIndexListEndOffset);
  }

  /* Write end signature. */
  fseek(f, 0L, SEEK_END);
  writeOne(f, sig);

  carefulClose(&f);
  freez(&bedString);
  freeHash(&lenHash);
  asObjectFree(&as);
  lineFileClose(&lf);
  bbiChromUsageFreeList(&usageList);
  popRHandlers();
  return r_outfile;
}
