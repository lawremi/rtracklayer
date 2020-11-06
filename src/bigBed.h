#ifndef BIG_BED_H
#define BIG_BED_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP BBDFile_seqlengths(SEXP r_filename);
SEXP BBDFile_fieldnames(SEXP r_filename);
SEXP BBDFile_query(SEXP r_filename, SEXP r_seqnames, SEXP r_ranges,
                   SEXP r_defaultindex, SEXP r_extraindex);
SEXP BBDFile_write(SEXP r_seqlengths, SEXP r_bedString, SEXP r_autosql,
                   SEXP r_indexfields, SEXP r_compress, SEXP r_outfile);

#endif
