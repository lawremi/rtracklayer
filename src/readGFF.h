#ifndef GFF_H
#define GFF_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP gff_colnames();
SEXP scan_gff(SEXP filexp, SEXP tags, SEXP filter);
SEXP load_gff(SEXP filexp, SEXP tags, SEXP filter,
	      SEXP ans_nrow, SEXP attrcol_fmt, SEXP pragmas,
	      SEXP colmap, SEXP raw_data);

#endif
