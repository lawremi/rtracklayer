/* basicBed.h contains the basic interface to Browser Extensible Data (bed) files and tables.
 * The idea behind bed is that the first three fields are defined and required.
 * A total of 15 fields are defined, and the file can contain any number of these.
 * In addition after any number of defined fields there can be custom fields that
 * are not defined in the bed spec.
 *
 * There's additional bed-related code in src/hg/inc/bed.h.  This module contains the
 * stuff that's independent of the database and other genomic structures. */

#ifndef BASICBED_H
#define BASICBED_H

#include "asParse.h"

struct bed
/* Browser extensible data */
    {
    struct bed *next;  /* Next in singly linked list. */
    char *chrom;	/* Human chromosome or FPC contig */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    char *name;	/* Name of item */

    /* The following items are not loaded by   the bedLoad routines. */
    int score; /* Score - 0-1000 */  /* Should be uint but there are still some ct users with neg values, .as DOES say uint */
    char strand[2];  /* + or -.  */
    unsigned thickStart; /* Start of where display should be thick (start codon for genes) */
    unsigned thickEnd;   /* End of where display should be thick (stop codon for genes) */
    unsigned itemRgb;    /* RGB 8 bits each */
    unsigned blockCount; /* Number of blocks. */
    int *blockSizes;     /* Comma separated list of block sizes.  */
    int *chromStarts;    /* Start positions inside chromosome.  Relative to chromStart*/


    int expCount;	/* Experiment count */
    int *expIds;		/* Comma separated list of Experiment ids */
    float *expScores;	/* Comma separated list of Experiment scores. */
    char *label;        /* Label to use on element if bigBed. */
    };

#define bedKnownFields 15	/* Maximum known fields in bed */

#define BB_MAX_CHROM_STRING 255  /* Maximum string length for chromosome length */

struct bed3
/* Browser extensible data - first three fields */
    {
    struct bed3 *next;  /* Next in singly linked list. */
    char *chrom;	/* Human chromosome or FPC contig */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    };


void bed3Free(struct bed3 **pBed);
/* Free up bed3 */


struct bed4
/* Browser extensible data - first four fields */
    {
    struct bed4 *next;  /* Next in singly linked list. */
    char *chrom;	/* Human chromosome or FPC contig */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    char *name;	/* Name of item */
    };


void bed4Free(struct bed4 **pBed);
/* Free up bed4 */


void bedFree(struct bed **pEl);
/* Free a single dynamically allocated bed such as created
 * with bedLoad(). */


#define bedTabOut(el,f) bedOutput(el,f,'\t','\n');
/* Print out bed as a line in a tab-separated file. */

#define bedCommaOut(el,f) bedOutput(el,f,',',',');
/* Print out bed as a comma separated list including final comma. */

/* --------------- End of AutoSQL generated code. --------------- */


struct bedLine
/* A line in a bed file with chromosome, start position parsed out. */
    {
    struct bedLine *next;	/* Next in list. */
    char *chrom;                /* Chromosome parsed out. */
    int chromStart;             /* Start position (still in rest of line). */
    char *line;                 /* Rest of line. */
    };

struct bedLine *bedLineNew(char *line);
/* Create a new bedLine based on tab-separated string s. */

void bedLineFree(struct bedLine **pBl);
/* Free up memory associated with bedLine. */

int bedLineCmp(const void *va, const void *vb);
/* Compare to sort based on chrom,chromStart. */

struct bed *bedLoad5(char **row);
/* Load first five fields of bed. */

struct bed *bedLoadN(char *row[], int wordCount);
/* Convert a row of strings to a bed. */

struct bed *bedLoadNAllChrom(char *fileName, int numFields, char* chrom);
/* Load bed entries from a tab-separated file that have the given chrom.
 * Dispose of this with bedFreeList(). */


void bedLoadAllReturnFieldCountAndRgb(char *fileName, struct bed **retList, int *retFieldCount, 
    boolean *retRgb);
/* Load bed of unknown size and return number of fields as well as list of bed items.
 * Ensures that all lines in bed file have same field count.  Also returns whether 
 * column 9 is being used as RGB or not. */


void bedOutFlexible(struct bed *el, int wordCount, FILE *f,
	char sep, char lastSep, boolean useItemRgb);
/* Write a bed of wordCount fields, optionally interpreting field nine as R,G,B values. */


int bedTotalBlockSize(struct bed *bed);
/* Return total size of all blocks. */


int bedBlockSizeInRange(struct bed *bed, int rangeStart, int rangeEnd);
/* Get size of all parts of all exons between rangeStart and rangeEnd. */


struct bed *cloneBed(struct bed *bed);
/* Make an all-newly-allocated copy of a single bed record. */


int bedParseRgb(char *itemRgb);
/*	parse a string: "r,g,b" into three unsigned char values
	returned as 24 bit number, or -1 for failure */


int bedSameStrandOverlap(struct bed *a, struct bed *b);
/* Return amount of block-level overlap on same strand between a and b */


struct rbTree *bedToRangeTree(struct bed *bed);
/* Convert bed into a range tree. */

void bedIntoRangeTree(struct bed *bed, struct rbTree *rangeTree);
/* Add all blocks in bed to range tree.  For beds without blocks,
 * add entire bed. */

int bedRangeTreeOverlap(struct bed *bed, struct rbTree *rangeTree);
/* Return number of bases bed overlaps with rangeTree. */

struct bed *bedThickOnly(struct bed *in);
/* Return a bed that only has the thick part. (Which is usually the CDS). */


char *bedAsDef(int bedFieldCount, int totalFieldCount);
/* Return an autoSql definition for a bed of given number of fields. 
 * Normally totalFieldCount is equal to bedFieldCount.  If there are extra
 * fields they are just given the names field16, field17, etc and type string. */


void loadAndValidateBedExt(char *row[], int bedFieldCount, int fieldCount, struct lineFile *lf, struct bed * bed, struct asObject *as, boolean isCt,  boolean allow1bpOverlap);
/* Convert a row of strings to a bed and validate the contents.  Abort with message if invalid data. Optionally validate bedPlus via asObject. Possibly allow one base overlap in exons */

int itemRgbColumn(char *column9);
/* Convert color specification to internal format. */
#endif /* BASICBED_H */
