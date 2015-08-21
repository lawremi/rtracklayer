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
	int i;

	x = strtod(val, &end_conversion);
	for (i = end_conversion - val; i < val_len; i++) {
		c = val[i];
		if (!isspace(c))
			return NA_REAL;
	}
	return x;
}

/* See http://www.sequenceontology.org/resources/gff3.html for the official
 * GFF3 specs. */
#define	SEQID_IDX 0
#define	SOURCE_IDX 1
#define	TYPE_IDX 2
#define	START_IDX 3
#define	END_IDX 4
#define	SCORE_IDX 5
#define	STRAND_IDX 6
#define	PHASE_IDX 7
#define	ATTRIBUTES_IDX 8

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

static SEXP alloc_ans(int nrow)
{
	SEXP ans, ans_names, ans_col, col_name;
	int col_idx;
	SEXPTYPE col_type;

	PROTECT(ans = NEW_LIST(9));
	PROTECT(ans_names = NEW_CHARACTER(9));
	for (col_idx = 0; col_idx < 9; col_idx++) {
		col_type = col_types[col_idx];
		PROTECT(ans_col = allocVector(col_type, nrow));
		SET_ELEMENT(ans, col_idx, ans_col);
		UNPROTECT(1);
		PROTECT(col_name = mkChar(col_names[col_idx]));
		SET_STRING_ELT(ans_names, col_idx, col_name);
		UNPROTECT(1);
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

static void collect_val(SEXP ans, int row_idx, int col_idx,
		const char *val, int val_len)
{
	SEXP ans_col;
	SEXPTYPE col_type;

	//printf("row_idx=%d col_idx=%d: ", row_idx, col_idx);
	ans_col = VECTOR_ELT(ans, col_idx);
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
		SEXP feature_types, SEXP ans, int *row_idx)
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
		if (ans != R_NilValue)
			collect_val(ans, *row_idx, col_idx,
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
	if (ans != R_NilValue) {
		val_len = delete_trailing_LF_or_CRLF(val, val_len);
		collect_val(ans, *row_idx, col_idx,
			    val, val_len);
	}
	(*row_idx)++;
	return NULL;
}

static const char *parse_GFF_file(SEXP filexp, SEXP feature_types,
		int *ans_nrow, SEXP ans)
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
			break;  /* stop parsing at first FASTA header */
		errmsg = parse_GFF_line(buf, lineno,
					feature_types, ans, &row_idx);
		if (errmsg != NULL)
			return errmsg;
	}
	*ans_nrow = row_idx;
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP GFFFile_read(SEXP filexp, SEXP feature_types)
{
	const char *errmsg;
	int ans_nrow;
	SEXP ans;

	errmsg = parse_GFF_file(filexp, feature_types,
				&ans_nrow, R_NilValue);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);
	PROTECT(ans = alloc_ans(ans_nrow));
	filexp_rewind(filexp);
	errmsg = parse_GFF_file(filexp, feature_types,
				&ans_nrow, ans);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);
	UNPROTECT(1);
	return ans;
}

