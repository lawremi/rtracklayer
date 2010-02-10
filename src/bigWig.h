#ifndef BIG_WIG_H
#define BIG_WIG_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP BWGSectionList_add(SEXP r_sections, SEXP r_seq, SEXP r_ranges,
                        SEXP r_score, SEXP r_format);
SEXP BWGSectionList_write(SEXP r_sections, SEXP r_seqlengths, SEXP r_compress,
                          SEXP r_file);
SEXP BWGSectionList_cleanup(SEXP r_sections);
SEXP BWGFile_query(SEXP r_filename, SEXP r_ranges, SEXP r_colnames);

#endif
