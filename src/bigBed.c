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
  struct asObject *as = bigBedAsOrDefault(file);
  struct asColumn *asCol = as->columnList;

  int fieldCount = file->fieldCount;
  int definedFieldCount = getDefinedFieldCount(as);
  bigBedFileClose(&file);

  SEXP defaultFields = PROTECT(allocVector(STRSXP, definedFieldCount));
  SEXP extraFields = PROTECT(allocVector(STRSXP, fieldCount - definedFieldCount));

  int extraIndex = 0;
  for (int i = 0; i < fieldCount; ++i) {
    if (i < definedFieldCount)
      SET_STRING_ELT(defaultFields, i, mkChar(asCol->name));
    else
      SET_STRING_ELT(extraFields, extraIndex++, mkChar(asCol->name));

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

static struct bigBedInterval *queryIntervals(struct bbiFile *file,
                                             SEXP r_seqnames,
                                             int *start, int *width,
                                             int n_ranges,
                                             SEXP n_qhits,
                                             struct lm *lm) {
    struct bigBedInterval *hits = NULL, *tail = NULL;
    for (int i = 0; i < n_ranges; ++i) {
        struct bigBedInterval *queryHits =
            bigBedIntervalQuery(file, (char *)CHAR(STRING_ELT(r_seqnames, i)),
                                start[i] - 1, start[i] - 1 + width[i], 0, lm);
        if (!hits) {
            hits = queryHits;
        } else {
            tail->next = queryHits;
        }
        if (queryHits)
            tail = slLastEl(queryHits);
        INTEGER(n_qhits)[i] = slCount(queryHits);
    }
    return hits;
}

static SEXPTYPE mapAsTypeToRType(enum asTypes ftype) {
    if (asTypesIsFloating(ftype))
        return REALSXP;

    switch (ftype) {
        case t_int:
        case t_short:
        case t_ushort:
        case t_byte:
        case t_ubyte:
        case t_uint:
        /* Assumes values fit in signed 32-bit range (typical genomic coords) */
            return INTSXP;

        /* 64-bit integer stored as double (no native 64-bit int in base R) */
        case t_off:
            return REALSXP;

        case t_char:
        case t_string:
        case t_lstring:
            return STRSXP;

        default:
            return STRSXP;
    }
}

static SEXP prepareFieldContainers(struct asObject *as,
                                   int fieldCount, int definedFieldCount,
                                   int n_hits, SEXP r_colnames,
                                   SEXPTYPE *fieldTypes,
                                   int *unprotectCount) {
    struct asColumn *asCol = as->columnList;
    SEXP fieldList = PROTECT(allocVector(VECSXP, fieldCount));
    SEXP fieldNames = PROTECT(allocVector(STRSXP, fieldCount));
    *unprotectCount += 2;

    struct hash *colHash = hashNew(0);
    int n_colnames = LENGTH(r_colnames);
    for (int k = 0; k < n_colnames; ++k) {
        hashAdd(colHash, (char *)CHAR(STRING_ELT(r_colnames, k)), NULL);
    }

    for (int j = 0; j < fieldCount; ++j, asCol = asCol->next) {
        char *colName = asCol->name;

        if (colName)
            SET_STRING_ELT(fieldNames, j, mkChar(colName));
        else
            SET_VECTOR_ELT(fieldNames, j, R_NilValue);

        bool selected = (colName && hashLookup(colHash, colName) != NULL);

        if (!selected) {
            SET_VECTOR_ELT(fieldList, j, R_NilValue);
            continue;
        }

        fieldTypes[j] = mapAsTypeToRType(asCol->lowType->type);
        SET_VECTOR_ELT(fieldList, j, allocVector(fieldTypes[j], n_hits));
        PROTECT(VECTOR_ELT(fieldList, j));
        ++(*unprotectCount);
    }
    setAttrib(fieldList, R_NamesSymbol, fieldNames);
    hashFree(&colHash);
    return fieldList;
}

static void fillResults(struct bigBedInterval *hits,
                        SEXP n_qhits, int n_ranges, int fieldCount,
                        int definedFieldCount, SEXPTYPE *fieldTypes,
                        SEXP fieldList, SEXP r_seqnames) {
    char startBuf[16], endBuf[16], *row[fieldCount];
    int i = 0, rangeIndex = 0, count = 0;

    for (struct bigBedInterval *bb = hits; bb; bb = bb->next, ++i, ++count) {
        if (rangeIndex < n_ranges && count == INTEGER(n_qhits)[rangeIndex]) {
            ++rangeIndex;
            count = 0;
        }

        bigBedIntervalToRow(bb, (char *)CHAR(STRING_ELT(r_seqnames, rangeIndex)),
                            startBuf, endBuf, row, fieldCount);

        for (int j = 0; j < fieldCount; ++j) {
            if (VECTOR_ELT(fieldList, j) == R_NilValue)
                continue;

            switch (fieldTypes[j]) {
                case REALSXP:
                    REAL(VECTOR_ELT(fieldList, j))[i] = sqlDouble(row[j]);
                    break;
                case INTSXP:
                    INTEGER(VECTOR_ELT(fieldList, j))[i] = sqlSigned(row[j]);
                    break;
                case STRSXP:
                    SET_STRING_ELT(VECTOR_ELT(fieldList, j), i, mkChar(row[j]));
                    break;
            }
        }
    }
}

static SEXP wrapResults(SEXP n_qhits, SEXP fieldList, int fieldCount,
                        int *unprotectCount) {
    int listSize = 2 + fieldCount;  // (n_qhits, fields)
    SEXP ans = PROTECT(allocVector(VECSXP, listSize));
    ++(*unprotectCount);

    int idx = 0;
    SET_VECTOR_ELT(ans, idx++, n_qhits);
    for (int j = 0; j < fieldCount; ++j) {
        SET_VECTOR_ELT(ans, idx++, VECTOR_ELT(fieldList, j));
    }

    SEXP ansNames = PROTECT(allocVector(STRSXP, listSize));
    ++(*unprotectCount);
    SET_STRING_ELT(ansNames, 0, mkChar("n_qhits"));

    SEXP fieldNames = getAttrib(fieldList, R_NamesSymbol);
    for (int j = 0; j < fieldCount; ++j) {
        SET_STRING_ELT(ansNames, j + 1, STRING_ELT(fieldNames, j));
    }

    setAttrib(ans, R_NamesSymbol, ansNames);
    return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_query(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges,
                   SEXP r_colnames) {
    pushRHandlers();

    int unprotectCount = 0;
    const char *fname = CHAR(asChar(r_filename));
    struct bbiFile *file = bigBedFileOpen((char *)fname);
    struct lm *lm = lmInit(0);

    int n_ranges = get_IRanges_length(r_ranges);
    int *start   = INTEGER(get_IRanges_start(r_ranges));
    int *width   = INTEGER(get_IRanges_width(r_ranges));

    SEXP n_qhits = PROTECT(allocVector(INTSXP, n_ranges));
    ++unprotectCount;

    struct bigBedInterval *hits = queryIntervals(file, r_seqnames, start, width,
                                                 n_ranges, n_qhits, lm);

    struct asObject *as   = bigBedAsOrDefault(file);
    int fieldCount        = file->fieldCount;
    int definedFieldCount = getDefinedFieldCount(as);
    int n_hits            = slCount(hits);

    bigBedFileClose(&file);

    SEXPTYPE *fieldTypes = (SEXPTYPE*)R_alloc(fieldCount, sizeof(SEXPTYPE));
    SEXP fieldList = prepareFieldContainers(as, fieldCount, definedFieldCount,
                                            n_hits, r_colnames, fieldTypes,
                                            &unprotectCount);

    fillResults(hits, n_qhits, n_ranges, fieldCount, definedFieldCount,
                fieldTypes, fieldList, r_seqnames);
    SEXP ans = wrapResults(n_qhits, fieldList, fieldCount, &unprotectCount);

    asObjectFree(&as);
    lmCleanup(&lm);
    UNPROTECT(unprotectCount);
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
#ifndef NDEBUG
  bits64 extraIndexListEndOffset = 0;
#endif
  if (extraIndexList != NULL) {
    extraIndexListOffset = ftell(f);
    int extraIndexSize = 16 + 4*1;   /* Fixed record size 16, plus 1 times field size of 4 */
    repeatCharOut(f, 0, extraIndexSize*extraIndexCount);
#ifndef NDEBUG
    extraIndexListEndOffset = ftell(f);
#endif
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
