#ifndef BIGBED_HELPER_H
#define BIGBED_HELPER_H

#include "ucsc/sig.h"
#include "ucsc/common.h"
#include "ucsc/sqlNum.h"
#include "ucsc/asParse.h"
#include "ucsc/obscure.h"
#include "ucsc/bbiFile.h"
#include "ucsc/zlibFace.h"
#include "ucsc/basicBed.h"
#include "ucsc/bPlusTree.h"
#include "ucsc/rangeTree.h"

#include "rtracklayer.h"

static int maxBedNameSize;

enum IFields
{
  i_name    = 3,  /* Index value of name field    */
  i_score   = 4,  /* Index value of score field   */
  i_strand  = 5,  /* Index value of strand field  */
  i_thick   = 7,  /* Index value of thick field   */
  i_itemRgb = 8,  /* Index value of itemRgb field */
  i_blocks  = 11, /* Index value of blocks field  */
};

int getDefinedFieldCount(struct asObject *as);
bool isPresent(int definedFieldCount, int index);
bool isSelected(SEXP r_selectedindex, int position);

void *bbNamedFileChunkVal(const void *va);
void bbNamedFileChunkKey(const void *va, char *keyBuf);
int bbNamedFileChunkCmpByName(const void *va, const void *vb);
struct rbTree *rangeTreeForBedChrom(struct lineFile *lf, char *chrom);
void bbExIndexMakerAllocChunkArrays(struct bbExIndexMaker *eim, int recordCount);
void bbExIndexMakerAddKeysFromRow(struct bbExIndexMaker *eim, char **row, int recordIx);
struct bbExIndexMaker *bbExIndexMakerNew(struct slName *extraIndexList, struct asObject *as);
void bbExIndexMakerAddOffsetSize(struct bbExIndexMaker *eim, bits64 offset, bits64 size,
                                 long startIx, long endIx);
struct bbiSummary *bedWriteReducedOnceReturnReducedTwice(struct bbiChromUsage *usageList,
    int fieldCount, struct lineFile *lf, bits32 initialReduction, bits32 initialReductionCount,
    int zoomIncrement, int blockSize, int itemsPerSlot, boolean doCompress,
    struct lm *lm, FILE *f, bits64 *retDataStart, bits64 *retIndexStart,
    struct bbiSummaryElement *totalSum);
void writeBlocks(struct bbiChromUsage *usageList, struct lineFile *lf, struct asObject *as,
    int itemsPerSlot, struct bbiBoundsArray *bounds,
    int sectionCount, boolean doCompress, FILE *f,
    int resTryCount, int resScales[], int resSizes[],
    struct bbExIndexMaker *eim,  int bedCount,
    bits16 fieldCount, int bedN, bits32 *retMaxBlockSize);

#endif
