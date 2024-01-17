/* Some stuff that you'll likely need in any program that works with
 * DNA.  Includes stuff for amino acids as well. 
 *
 * Assumes that DNA is stored as a character.
 * The DNA it generates will include the bases 
 * as lowercase tcag.  It will generally accept
 * uppercase as well, and also 'n' or 'N' or '-'
 * for unknown bases. 
 *
 * Amino acids are stored as single character upper case. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */


#ifndef DNAUTIL_H
#define DNAUTIL_H

void dnaUtilOpen(); /* Good idea to call this before using any arrays
		     * here.  */

/* Numerical values for bases. */
#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4   /* Used in 1/2 byte representation. */

typedef char DNA;
typedef char AA;
typedef char BIOPOL;	/* Biological polymer. */

/* A little array to help us decide if a character is a 
 * nucleotide, and if so convert it to lower case. 
 * Contains zeroes for characters that aren't used
 * in DNA sequence. */
extern DNA ntChars[256];
extern AA aaChars[256];

/* An array that converts alphabetical DNA representation
 * to numerical one: X_BASE_VAL as above.  For charaters
 * other than [atgcATGC], has -1. */
extern int ntVal[256];
extern int aaVal[256];
extern int ntValLower[256];	/* NT values only for lower case. */
extern int ntValUpper[256];	/* NT values only for upper case. */

/* Like ntVal, but with T_BASE_VAL in place of -1 for nonexistent nucleotides. */
extern int ntValNoN[256];     

/* Like ntVal but with N_BASE_VAL in place of -1 for 'n', 'x', '-', etc. */
extern int ntVal5[256];

/* Inverse array - takes X_BASE_VAL int to a DNA char
 * value. */
extern DNA valToNt[];

/* Similar array that doesn't convert to lower case. */
extern DNA ntMixedCaseChars[256];

/* Another array to help us do complement of DNA  */
extern DNA ntCompTable[256];

/* Arrays to convert between lower case indicating repeat masking, and
 * a 1/2 byte representation where the 4th bit indicates if the characeter
 * is masked. Uses N_BASE_VAL for `n', `x', etc.
*/
#define MASKED_BASE_BIT 8
extern int ntValMasked[256];
extern DNA valToNtMasked[256];

/*Complement DNA (not reverse)*/
void complement(DNA *dna, long length);

/* Reverse complement DNA. */
void reverseComplement(DNA *dna, long length);

enum dnaCase {dnaUpper,dnaLower,dnaMixed,};
/* DNA upper, lower, or mixed case? */

typedef char Codon; /* Our codon type. */

/* Return single letter code (upper case) for protein.
 * Returns X for bad input, 0 for stop codon.
 * The "Standard" Code */
AA lookupCodon(DNA *dna); 


/* Returns one letter code for protein, 
 * 0 for stop codon or X for bad input,
 * Vertebrate Mitochondrial Code */
AA lookupMitoCodon(DNA *dna);


extern char *aaAbbr(int i);
/* return pointer to AA abbrevation */

extern char aaLetter(int i);
/* return AA letter */

UBYTE packDna4(DNA *in);
/* Pack 4 bases into a UBYTE */

void unalignedUnpackDna(bits32 *tiles, int start, int size, DNA *unpacked);
/* Unpack into out, even though not starting/stopping on tile 
 * boundaries. */

int intronOrientationMinSize(DNA *iStart, DNA *iEnd, int minIntronSize);
/* Given a gap in genome from iStart to iEnd, return 
 * Return 1 for GT/AG intron between left and right, -1 for CT/AC, 0 for no
 * intron. */


int  dnaOrAaScoreMatch(char *a, char *b, int size, int matchScore, int mismatchScore, 
	char ignore);
/* Compare two sequences (without inserts or deletions) and score. */


boolean isDna(char *poly, int size);
/* Return TRUE if letters in poly are at least 90% ACGTU */


#endif /* DNAUTIL_H */
