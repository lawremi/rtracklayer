/* dnaSeq.c - stuff to manage DNA sequences. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "dnaseq.h"
#include "bits.h"
#include "hash.h"
#include "obscure.h"



struct dnaSeq *newDnaSeq(DNA *dna, int size, char *name)
/* Create a new DNA seq. */
{
struct dnaSeq *seq;

seq = needMem(sizeof(*seq));
if (name != NULL)
    seq->name = cloneString(name);
seq->dna = dna;
seq->size = size;
seq->mask = NULL;
return seq;
}

void freeDnaSeq(struct dnaSeq **pSeq)
/* Free up DNA seq. (And unlink underlying resource node.) */
{
struct dnaSeq *seq = *pSeq;
if (seq == NULL)
    return;
freeMem(seq->name);
freeMem(seq->dna);
bitFree(&seq->mask);
freez(pSeq);
}


aaSeq *translateSeqN(struct dnaSeq *inSeq, unsigned offset, unsigned inSize, boolean stop)
/* Return a translated sequence.  Offset is position of first base to
 * translate. If size is 0 then use length of inSeq. */
{
aaSeq *seq;
DNA *dna = inSeq->dna;
AA *pep, aa;
int i, lastCodon;
int actualSize = 0;

assert(offset <= inSeq->size);
if ((inSize == 0) || (inSize > (inSeq->size - offset)))
    inSize = inSeq->size - offset;
lastCodon = offset + inSize - 3;

AllocVar(seq);
seq->dna = pep = needLargeMem(inSize/3+1);
for (i=offset; i <= lastCodon; i += 3)
    {
    aa = lookupCodon(dna+i);
    if (aa == 0)
	{
        if (stop)
	    break;
	else
	    aa = 'Z';
	}
    *pep++ = aa;
    ++actualSize;
    }
*pep = 0;
assert(actualSize <= inSize/3+1);
seq->size = actualSize;
seq->name = cloneString(inSeq->name);
return seq;
}
