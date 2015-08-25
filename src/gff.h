#ifndef GFF_H
#define GFF_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP gff_colnames();
SEXP gff_read(SEXP filexp, SEXP colmap, SEXP tags, SEXP filter, SEXP raw_data);

#endif
