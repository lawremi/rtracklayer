#ifndef TWO_BIT_H
#define TWO_BIT_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP DNAString_to_twoBit(SEXP r_dna, SEXP r_mask, SEXP r_seqname);
SEXP TwoBits_write(SEXP r_twoBits, SEXP r_filename);
SEXP TwoBitFile_seqlengths(SEXP r_filename);
SEXP TwoBitFile_read(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges,
                     SEXP lkup);

#endif
