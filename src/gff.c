#include "XVector_interface.h"
#include "S4Vectors_interface.h"

#include <ctype.h>   /* for isdigit() and isspace() */
#include <stdlib.h>  /* for strtod() */

/*
 * Turn string pointed by 'val' into an int. The string has no terminating
 * null byte ('\0') and must have the following format:
 *     ^[[:space:]]*[+-]?[[:digit:]]+[[:space:]]*$
 * Return NA_INTEGER if the string is malformed or if it represents a integer
 * value that cannot be represented by an int (int overflow).
 */
#define	LEADING_SPACE 0
#define	NUMBER 1
#define	TRAILING_SPACE 2
static int as_int(const char *val, int val_len)
{
	int n, ndigit, sign, status, i;
	char c;

	n = ndigit = 0;
	sign = 1;
	status = LEADING_SPACE;
	for (i = 0; i < val_len; i++) {
		c = val[i];
		//printf("%c", c);
		if (isdigit(c)) {
			if (status == TRAILING_SPACE)
				return NA_INTEGER;  /* malformed string */
			status = NUMBER;
			ndigit++;
			n = safe_int_mult(n, 10);
			n = safe_int_add(n, c - '0');
			if (n == NA_INTEGER)
				return NA_INTEGER;  /* int overflow */
			continue;
		}
		if (c == '+' || c == '-') {
			if (status != LEADING_SPACE)
				return NA_INTEGER;  /* malformed string */
			status = NUMBER;
			if (c == '-')
				sign = -1;
			continue;
		}
		if (!isspace(c))
			return NA_INTEGER;  /* malformed string */
		if (status == NUMBER) {
			if (ndigit == 0)
				return NA_INTEGER;  /* malformed string */
			status = TRAILING_SPACE;
		}
	}
	if (ndigit == 0)
		return NA_INTEGER;  /* malformed string */
	if (sign == -1)
		n = -n;
	return n;
}

static double as_double(const char *val, int val_len)
{
	double x;
	char *end_conversion, c;
	int end_offset, i;

	x = strtod(val, &end_conversion);
	end_offset = end_conversion - val;
	if (end_offset == 0)
		return NA_REAL;
	for (i = end_offset; i < val_len; i++) {
		c = val[i];
		if (!isspace(c))
			return NA_REAL;
	}
	return x;
}

/* See http://www.sequenceontology.org/resources/gff3.html for the official
 * GFF3 specs. */

static const char *col_names[] = {
	"seqid",
	"source",
	"type",
	"start",
	"end",
	"score",
	"strand",
	"phase",
	"attributes"
};

static const SEXPTYPE col_types[] = {
	STRSXP,   /* seqid */
	STRSXP,   /* source */
	STRSXP,   /* type */
	INTSXP,   /* start */
	INTSXP,   /* end */
	REALSXP,  /* score */
	STRSXP,   /* strand */
	INTSXP,   /* phase */
	STRSXP    /* attributes */
};

#define	GFF_NCOL (sizeof(col_names) / sizeof(char *))

#define	SEQID_IDX 0
#define	SOURCE_IDX 1
#define	TYPE_IDX 2
#define	START_IDX 3
#define	END_IDX 4
#define	SCORE_IDX 5
#define	STRAND_IDX 6
#define	PHASE_IDX 7
#define	ATTRIBUTES_IDX 8

static int prepare_col_map(int *col_map, SEXP cols)
{
	int ncol, col_idx;

	if (!(cols == R_NilValue ||
	      (IS_LOGICAL(cols) && LENGTH(cols) == GFF_NCOL)))
		error("'cols' must be NULL or a logical vector of length %d",
		      GFF_NCOL);
	ncol = 0;
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		if (cols == R_NilValue || LOGICAL(cols)[col_idx]) {
			col_map[col_idx] = ncol++;
		} else {
			col_map[col_idx] = NA_INTEGER;
		}
		//printf("col_map[%d] = %d\n", col_idx, col_map[col_idx]);
	}
	return ncol;
}

static SEXP alloc_ans(int nrow, int ncol, const int *col_map)
{
	SEXP ans, ans_names, ans_col, col_name;
	int col_idx, j;
	SEXPTYPE col_type;

	PROTECT(ans = NEW_LIST(ncol));
	PROTECT(ans_names = NEW_CHARACTER(ncol));
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		j = col_map[col_idx];
		if (j == NA_INTEGER)
			continue;
		col_type = col_types[col_idx];
		PROTECT(ans_col = allocVector(col_type, nrow));
		SET_ELEMENT(ans, j, ans_col);
		UNPROTECT(1);
		PROTECT(col_name = mkChar(col_names[col_idx]));
		SET_STRING_ELT(ans_names, j, col_name);
		UNPROTECT(1);
		j++;
	}
	SET_NAMES(ans, ans_names);
	/* list_as_data_frame() performs IN-PLACE coercion */
	list_as_data_frame(ans, nrow);
	UNPROTECT(2);
	return ans;
}

static void collect_string(SEXP ans_col, int row_idx,
		const char *val, int val_len)
{
	SEXP tmp;

	if (val_len == 1 && val[0] == '.') {
		SET_STRING_ELT(ans_col, row_idx, NA_STRING);
	} else {
		PROTECT(tmp = mkCharLen(val, val_len));
		SET_STRING_ELT(ans_col, row_idx, tmp);
		UNPROTECT(1);
	}
	return;
}

static void collect_int(SEXP ans_col, int row_idx,
		const char *val, int val_len)
{
	INTEGER(ans_col)[row_idx] = as_int(val, val_len);
	return;
}

static void collect_double(SEXP ans_col, int row_idx,
		const char *val, int val_len)
{
	REAL(ans_col)[row_idx] = as_double(val, val_len);
	return;
}

static void collect_val(SEXP ans, int row_idx, int col_idx, const int *col_map,
		const char *val, int val_len)
{
	SEXP ans_col;
	SEXPTYPE col_type;

	//printf("row_idx=%d col_idx=%d: ", row_idx, col_idx);
	ans_col = VECTOR_ELT(ans, col_map[col_idx]);
	col_type = col_types[col_idx];
	switch (col_type) {
	    case STRSXP:
		if (col_idx == STRAND_IDX
		 && val_len == 1 && (val[0] == '.' || val[0] == '?'))
		{
			val = "*";
			val_len = 1;
		}
		collect_string(ans_col, row_idx, val, val_len);
	    break;
	    case INTSXP:
		collect_int(ans_col, row_idx, val, val_len);
	    break;
	    case REALSXP:
		collect_double(ans_col, row_idx, val, val_len);
	    break;
	}
	//printf("\n");
	return;
}

#define IOBUF_SIZE 20002
static char errmsg_buf[200];

static const char *parse_GFF_line(const char *line, int lineno,
		SEXP feature_types,
		SEXP ans, int *row_idx, const int *col_map)
{
	int col_idx, i, val_len;
	const char *val;
	char c;

	col_idx = i = 0;
	val = line;
	val_len = 0;
	while ((c = line[i++])) {
		if (c != '\t') {
			val_len++;
			continue;
		}
		if (col_idx >= 8) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "line %d has more than 8 tabs", lineno);
			return errmsg_buf;
		}
		if (ans != R_NilValue && col_map[col_idx] != NA_INTEGER)
			collect_val(ans, *row_idx, col_idx, col_map,
				    val, val_len);
		col_idx++;
		val = line + i;
		val_len = 0;
	}
	if (col_idx < 7) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "line %d has less than 7 tabs", lineno);
		return errmsg_buf;
	}
	if (ans != R_NilValue && col_map[col_idx] != NA_INTEGER) {
		val_len = delete_trailing_LF_or_CRLF(val, val_len);
		collect_val(ans, *row_idx, col_idx, col_map,
			    val, val_len);
	}
	(*row_idx)++;
	return NULL;
}

static const char *parse_GFF_file(SEXP filexp, SEXP feature_types,
		int *ans_nrow, SEXP ans, const int *col_map)
{
	int row_idx, lineno, ret_code, EOL_in_buf;
	char buf[IOBUF_SIZE], c;
	const char *errmsg;

	row_idx = 0;
	for (lineno = 1;
	     (ret_code = filexp_gets(filexp, buf, IOBUF_SIZE, &EOL_in_buf));
	     lineno += EOL_in_buf)
	{
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (!EOL_in_buf) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, "
				 "line is too long", lineno);
			return errmsg_buf;
		}
		c = buf[0];
		if (c == '\n' || (c == '\r' && buf[1] == '\n'))
			continue;  /* skip empty line */
		if (c == '#')
			continue;  /* skip comment */
		if (c == '>')
			break;     /* stop parsing at first FASTA header */
		errmsg = parse_GFF_line(buf, lineno, feature_types,
					ans, &row_idx, col_map);
		if (errmsg != NULL)
			return errmsg;
	}
	*ans_nrow = row_idx;
	return NULL;
}

/*
 * --- .Call ENTRY POINT ---
 * Args:
 *   filexp:        A "File External Pointer" (see src/io_utils.c in the
 *                  XVector package).
 *   cols:          NULL or a logical vector of length 9 indicating the
 *                  columns to collect. NULL means collect all columns.
 *   feature_types: NULL or a character vector of feature types. Only rows
 *                  with that type are collected. NULL means collect all rows.
 */
SEXP GFFFile_read(SEXP filexp, SEXP cols, SEXP feature_types)
{
	int col_map[GFF_NCOL], ans_nrow, ans_ncol;
	const char *errmsg;
	SEXP ans;

	ans_ncol = prepare_col_map(col_map, cols);
	if (!(feature_types == R_NilValue || IS_CHARACTER(feature_types)))
		error("'feature_types' must be NULL or a character vector");

	/* 1st pass */
	errmsg = parse_GFF_file(filexp, feature_types,
				&ans_nrow, R_NilValue, col_map);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);

	/* 2nd pass */
	PROTECT(ans = alloc_ans(ans_nrow, ans_ncol, col_map));
	filexp_rewind(filexp);
	errmsg = parse_GFF_file(filexp, feature_types,
				&ans_nrow, ans, col_map);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);
	UNPROTECT(1);
	return ans;
}

