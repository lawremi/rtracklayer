/* dnaSeq - stuff to manage DNA sequences. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef DNASEQ_H
#define DNASEQ_H

#ifndef DNAUTIL_H
#include "dnautil.h"
#endif

#ifndef BITS_H
#include "bits.h"
#endif

struct dnaSeq
/* A dna sequence in one-character per base format. */
    {
    struct dnaSeq *next;  /* Next in list. */
    char *name;           /* Name of sequence. */
    DNA *dna;             /* Sequence base by base. */
    int size;             /* Size of sequence. */
    Bits* mask;           /* Repeat mask (optional) */
    };

typedef struct dnaSeq bioSeq;	/* Preferred use if either DNA or protein. */
typedef struct dnaSeq aaSeq;	/* Preferred use if protein. */

struct dnaSeq *newDnaSeq(DNA *dna, int size, char *name);
/* Create a new DNA seq. */

void freeDnaSeq(struct dnaSeq **pSeq);
/* Free up DNA seq.  */
#define dnaSeqFree freeDnaSeq

aaSeq *translateSeqN(struct dnaSeq *inSeq, unsigned offset, unsigned size, boolean stop);
/* Return a translated sequence.  Offset is position of first base to
 * translate. If size is 0 then use length of inSeq. */


#endif /* DNASEQ_H */

