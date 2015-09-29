#ifndef GFF_H
#define GFF_H

#include "rtracklayer.h"

/* The .Call entry points */

SEXP gff_colnames(SEXP GFF1);
SEXP read_gff_pragmas(SEXP filexp);
SEXP scan_gff(SEXP filexp, SEXP attrcol_fmt, SEXP tags,
	      SEXP filter, SEXP nrows);
SEXP load_gff(SEXP filexp, SEXP attrcol_fmt, SEXP tags,
	      SEXP filter, SEXP nrows,
	      SEXP pragmas, SEXP colmap, SEXP raw_data);

#endif
