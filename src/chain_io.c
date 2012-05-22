#include <S.h>
#include "ucsc/common.h"
#include "ucsc/hash.h"

#include "rtracklayer.h"

#define LINEBUF_SIZE 20001

/* hash these chain blocks by target and query name */
typedef struct _ChainBlock {
  char *name;
  RangeAE ranges; /* to become an IRanges */
  IntAE offset; /* starts in the other sequence */
  /* rle of spaces and scores */
  IntAE length, score;
  CharAE rev; /* use CharAE until we have a bitset */
  CharAEAE space;
} ChainBlock;

#define HEADER_SIZE 11
#define DATA_SIZE 3

/* returns an array of ChainBlock pointers */
ChainBlock **read_chain_file(FILE *stream, const char *exclude, int *nblocks) {
  /* fgets() a line
     if first or after blank, parse header
       get names, use hash table to get existing block, ow create one
       add score, get offsets for starts
     ow parse a record for width, tstart, qstart
  */
  char linebuf[LINEBUF_SIZE];
  char *header[HEADER_SIZE];
  char *data[DATA_SIZE];
  int tstart, qstart;
  Rboolean new_block = TRUE, excluded = FALSE, trc, qrc;
  ChainBlock *block, **result;
  struct hash *hash = hashNew(6);
  struct hashEl *hash_elements;
  int line = 0, i = 0, header_line;
  while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
    line++;
    if (strlen(linebuf) == LINEBUF_SIZE - 1) {
      error("line %d is too long", line);
    }
    if (excluded) {
      eraseWhiteSpace(linebuf);
      if (!strlen(linebuf)) {
        excluded = FALSE;
        new_block = TRUE;
      }
    } else if (new_block) { /* have a header */
      void *value;
      int matches = chopByChar(linebuf, ' ', header, HEADER_SIZE);
      if (matches < HEADER_SIZE)
        error("expected %d elements in header, got %d, on line %d", HEADER_SIZE,
              matches, line);
      new_block = FALSE;
      if (exclude && (strstr(header[2], exclude) || strstr(header[7], exclude)))
        {
          //Rprintf("excluding: [%s -> %s]\n", header[2], header[7]);
          excluded = TRUE;
          continue;
        }
      value = hashFindVal(hash, header[2]);
      if (!value) { /* new block */
        int name_size = strlen(header[2])+1;
        block = Salloc(1, ChainBlock);
        hashAdd(hash, header[2], block);
        block->name = Salloc(name_size, char);
        memcpy(block->name, header[2], name_size);
        block->ranges = new_RangeAE(0, 0);
        block->offset = new_IntAE(0, 0, 0);
        block->length = new_IntAE(0, 0, 0);
        block->score = new_IntAE(0, 0, 0);
        block->rev = new_CharAE(0);
        block->space = new_CharAEAE(0, 0);
      } else block = value;
      IntAE_insert_at(&block->score, IntAE_get_nelt(&block->score),
                      atoi(header[1]));
      append_string_to_CharAEAE(&block->space, header[7]);
      header_line = line;
      trc = strcmp("+", header[4]);
      qrc = strcmp("+", header[9]);
      CharAE_insert_at(&block->rev, CharAE_get_nelt(&block->rev),
                       trc != qrc);
      tstart = atoi(header[5]) + 1; /* 0-based -> 1-based */
      if (trc)
        tstart = atoi(header[3]) - tstart + 2; /* start one too high */
      qstart = atoi(header[10]) + 1;
      if (qrc)
        qstart = atoi(header[8]) - qstart + 2;
    } else {
      int matches = chopByChar(linebuf, '\t', data, DATA_SIZE), width;
      if (matches != 1 && matches != 3)
        error("expecting 1 or 3 elements on line %d, got %d", line, matches);
      width = atoi(data[0]);
      tstart -= (trc ? width : 0);
      qstart -= (qrc ? width : 0);
      RangeAE_insert_at(&block->ranges, RangeAE_get_nelt(&block->ranges),
                        tstart, width);
      IntAE_insert_at(&block->offset, IntAE_get_nelt(&block->offset),
                      tstart - qstart);
      if (matches == 3) { /* normal line */
        int dt = atoi(data[1]), dq = atoi(data[2]);
        int tchange, qchange;
        if (trc) /* width already subtracted above */
          tchange = -dt;
        else tchange = width + dt;
        tstart += tchange;
        if (qrc)
          qchange = -dq;
        else qchange = width + dq;
        qstart += qchange;
      } else {
        new_block = TRUE;
        IntAE_insert_at(&block->length, IntAE_get_nelt(&block->length),
                        line-header_line);
        //Rprintf("end of %s block, line: %d\n", block->name, line);
        fgets(linebuf, LINEBUF_SIZE, stream); /* skip empty line */
        line++;
      }
    }
  }
  result = Salloc(hashNumEntries(hash), ChainBlock *);
  hash_elements = hashElListHash(hash);
  for (struct hashEl *h = hash_elements; h; h = h->next, i++) {
    result[i] = h->val;
  }
  *nblocks = i;
  hashElFreeList(&hash_elements);
  hashFree(&hash);
  return result;
}

/* parse the chain file format from UCSC for representing alignments */
SEXP readChain(SEXP r_path, SEXP r_exclude) {
  const char *path, *exclude;
  FILE *stream;
  SEXP ans, ans_names, ans_listData, chain_class, chainBlock_class;
  ChainBlock **chains;
  int i, nblocks;

  path = translateChar(STRING_ELT(r_path, 0));
  if ((stream = fopen(path, "r")) == NULL)
    error("cannot open file '%s'", path);
  exclude = r_exclude == R_NilValue ? NULL : CHAR(STRING_ELT(r_exclude, 0));
  chains = read_chain_file(stream, exclude, &nblocks);

  PROTECT(chain_class = MAKE_CLASS("Chain"));
  PROTECT(chainBlock_class = MAKE_CLASS("ChainBlock"));
  
  PROTECT(ans = NEW_OBJECT(chain_class));
  ans_listData = allocVector(VECSXP, nblocks);
  SET_SLOT(ans, install("listData"), ans_listData);
  ans_names = allocVector(STRSXP, nblocks);
  SET_NAMES(ans_listData, ans_names);
  for (i = 0; i < nblocks; i++) {
    SEXP block;
    block = NEW_OBJECT(chainBlock_class);
    SET_VECTOR_ELT(ans_listData, i, block);
    SET_SLOT(block, install("ranges"),
		new_IRanges_from_RangeAE("IRanges", &chains[i]->ranges));
    SET_SLOT(block, install("offset"),
		new_INTEGER_from_IntAE(&chains[i]->offset));
    SET_SLOT(block, install("length"),
		new_INTEGER_from_IntAE(&chains[i]->length));
    SET_SLOT(block, install("score"),
		new_INTEGER_from_IntAE(&chains[i]->score));
    SET_SLOT(block, install("space"),
		new_CHARACTER_from_CharAEAE(&chains[i]->space));
    SET_SLOT(block, install("reversed"),
		new_LOGICAL_from_CharAE(&chains[i]->rev));
    SET_STRING_ELT(ans_names, i, mkChar(chains[i]->name));
  }

  UNPROTECT(3);

  return ans;
}
