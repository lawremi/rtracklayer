#ifndef TWO_BIT_H
#define TWO_BIT_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP DNAString_to_twoBit(SEXP r_dna, SEXP r_mask, SEXP r_seqname);

#endif
