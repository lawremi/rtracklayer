/* Stuff that should be in ucsc/bwgInternal.h, but wasn't, so
   rtracklayer put it here. */

struct bwgBedGraphItem
/* An bedGraph-type item in a bwgSection. */
{
  struct bwgBedGraphItem *next;	/* Next in list. */
  bits32 start,end;		/* Range of chromosome covered. */
  float val;			/* Value. */
};

struct bwgVariableStepItem
/* An variableStep type item in a bwgSection. */
{
  struct bwgVariableStepItem *next;	/* Next in list. */
  bits32 start;		/* Start position in chromosome. */
  float val;			/* Value. */
};

struct bwgVariableStepPacked
/* An variableStep type item in a bwgSection. */
{
  bits32 start;		/* Start position in chromosome. */
  float val;			/* Value. */
};

struct bwgFixedStepItem
/* An fixedStep type item in a bwgSection. */
{
  struct bwgFixedStepItem *next;	/* Next in list. */
  float val;			/* Value. */
};

struct bwgFixedStepPacked
/* An fixedStep type item in a bwgSection. */
{
  float val;			/* Value. */
};

union bwgItem
/* Union of item pointers for all possible section types. */
{
  struct bwgBedGraphItem *bedGraphList;		/* A linked list */
  struct bwgFixedStepPacked *fixedStepPacked;		/* An array */
  struct bwgVariableStepPacked *variableStepPacked;	/* An array */
  /* No packed format for bedGraph... */
};

struct bwgSection
/* A section of a bigWig file - all on same chrom.  This is a somewhat fat data
 * structure used by the bigWig creation code.  See also bwgSection for the
 * structure returned by the bigWig reading code. */
{
  struct bwgSection *next;		/* Next in list. */
  char *chrom;			/* Chromosome name. */
  bits32 start,end;			/* Range of chromosome covered. */
  enum bwgSectionType type;
  union bwgItem items;		/* List/array of items in this section. */
  bits32 itemStep;			/* Step within item if applicable. */
  bits32 itemSpan;			/* Item span if applicable. */
  bits16 itemCount;			/* Number of items in section. */
  bits32 chromId;			/* Unique small integer value for chromosome. */
  bits64 fileOffset;			/* Offset of section in file. */
};

void bwgCreate(struct bwgSection *sectionList, struct hash *chromSizeHash, 
               int blockSize, int itemsPerSlot, boolean doCompress,
               char *fileName);

struct bbiFile *bigWigFileOpen(char *fileName);

struct bbiInterval *bigWigIntervalQuery(struct bbiFile *bwf, char *chrom,
                                        bits32 start, bits32 end,
                                        struct lm *lm);
