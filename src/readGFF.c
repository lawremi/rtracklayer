#include "XVector_interface.h"
#include "S4Vectors_interface.h"

#include <ctype.h>   /* for isdigit() and isspace() */
#include <stdlib.h>  /* for strtod() */
#include <string.h>  /* for memcpy() and memcmp() */

/*
#include <time.h>
static clock_t clock0;
static void init_clock(const char *msg)
{
	printf("%s", msg);
	clock0 = clock();
}
static void print_elapsed_time()
{
	printf("%8.6f s\n", ((double) clock() - clock0) / CLOCKS_PER_SEC);
}
*/

/*
 * Turn string pointed by 'val' into an int. The string has no terminating
 * null byte ('\0') and must have the following format:
 *     ^[[:space:]]*[+-]?[[:digit:]]+[[:space:]]*$
 * Return NA_INTEGER if the string is malformed or if it represents a integer
 * value that cannot be represented by an int (int overflow).
 * TODO: Maybe implement this on top of strtol(). Would be much simpler but
 * would it be equivalent? Also would it be as fast? See how as_double() below
 * is implemented on top of strtod().
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

#define	GFF_NCOL ((int) (sizeof(col_names) / sizeof(char *)))

/* --- .Call ENTRY POINT --- */
SEXP gff_colnames()
{
	SEXP ans, ans_elt;
	int col_idx;

	PROTECT(ans = NEW_CHARACTER(GFF_NCOL));
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		PROTECT(ans_elt = mkChar(col_names[col_idx]));
		SET_STRING_ELT(ans, col_idx, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

#define	SEQID_IDX 0
#define	SOURCE_IDX 1
#define	TYPE_IDX 2
#define	START_IDX 3
#define	END_IDX 4
#define	SCORE_IDX 5
#define	STRAND_IDX 6
#define	PHASE_IDX 7
#define	ATTRIBUTES_IDX 8

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

static int prepare_colmap0(int *colmap0, SEXP colmap)
{
	int ans_ncol0, col_idx, j;

	ans_ncol0 = 0;
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		j = INTEGER(colmap)[col_idx];
		if (j != NA_INTEGER) {
			if (j > ans_ncol0)
				ans_ncol0 = j;
			j--;
		}
		colmap0[col_idx] = j;
	}
	return ans_ncol0;
}

static SEXP alloc_ans(int ans_nrow, int ans_ncol0, const int *colmap0,
		SEXP tags, SEXP pragmas, SEXP raw_data)
{
	int ans_ntag, ans_ncol, is_raw, col_idx, j, i;
	SEXP ans, ans_attr, ans_names, ans_col, col_name, tags_elt;
	SEXPTYPE col_type;

	ans_ntag = LENGTH(tags);
	ans_ncol = ans_ncol0 + ans_ntag;
	is_raw = LOGICAL(raw_data)[0];

	PROTECT(ans = NEW_LIST(ans_ncol));
	PROTECT(ans_names = NEW_CHARACTER(ans_ncol));

	/* Alloc standard GFF columns. */
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		j = colmap0[col_idx];
		if (j == NA_INTEGER)
			continue;
		col_type = is_raw ? STRSXP : col_types[col_idx];
		PROTECT(ans_col = allocVector(col_type, ans_nrow));
		SET_ELEMENT(ans, j, ans_col);
		UNPROTECT(1);
		PROTECT(col_name = mkChar(col_names[col_idx]));
		SET_STRING_ELT(ans_names, j, col_name);
		UNPROTECT(1);
		j++;
	}

	/* Alloc tag columns. */
	for (j = ans_ncol0; j < ans_ncol; j++) {
		PROTECT(ans_col = NEW_CHARACTER(ans_nrow));
		for (i = 0; i < ans_nrow; i++)
			SET_STRING_ELT(ans_col, i, NA_STRING);
		SET_ELEMENT(ans, j, ans_col);
		UNPROTECT(1);
		tags_elt = STRING_ELT(tags, j - ans_ncol0);
		PROTECT(col_name = duplicate(tags_elt));
		SET_STRING_ELT(ans_names, j, col_name);
		UNPROTECT(1);
	}

	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* list_as_data_frame() performs IN-PLACE coercion */
	list_as_data_frame(ans, ans_nrow);

	/* Set additional attributes. */
	PROTECT(ans_attr = ScalarInteger(ans_ncol0));
	SET_ATTR(ans, install("ncol0"), ans_attr);
	UNPROTECT(1);
	PROTECT(ans_attr = ScalarInteger(ans_ntag));
	SET_ATTR(ans, install("ntag"), ans_attr);
	UNPROTECT(1);
	PROTECT(ans_attr = duplicate(raw_data));
	SET_ATTR(ans, install("raw_data"), raw_data);
	UNPROTECT(1);
	PROTECT(ans_attr = duplicate(pragmas));
	SET_ATTR(ans, install("pragmas"), ans_attr);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

static void load_string(const char *data, int data_len,
		SEXP ans_col, int row_idx)
{
	SEXP tmp;

	PROTECT(tmp = mkCharLen(data, data_len));
	SET_STRING_ELT(ans_col, row_idx, tmp);
	UNPROTECT(1);
	return;
}

static void load_int(const char *data, int data_len,
		SEXP ans_col, int row_idx)
{
	INTEGER(ans_col)[row_idx] = as_int(data, data_len);
	return;
}

static void load_double(const char *data, int data_len,
		SEXP ans_col, int row_idx)
{
	REAL(ans_col)[row_idx] = as_double(data, data_len);
	return;
}

static void load_data(const char *data, int data_len,
		SEXP ans, int row_idx, int col_idx, const int *colmap0)
{
	SEXP ans_col;
	int is_raw;
	SEXPTYPE col_type;

	ans_col = VECTOR_ELT(ans, colmap0[col_idx]);
	is_raw = LOGICAL(GET_ATTR(ans, install("raw_data")))[0];
	if (is_raw) {
		load_string(data, data_len, ans_col, row_idx);
		return;
	}
	col_type = col_types[col_idx];
	switch (col_type) {
	    case STRSXP:
		if (data_len == 1) {
			if (col_idx == STRAND_IDX
			 && (data[0] == '.' || data[0] == '?'))
			{
				data = "*";
				data_len = 1;
			} else if (data[0] == '.') {
				SET_STRING_ELT(ans_col, row_idx, NA_STRING);
				break;
			}
		}
		load_string(data, data_len, ans_col, row_idx);
	    break;
	    case INTSXP:
		load_int(data, data_len, ans_col, row_idx);
	    break;
	    case REALSXP:
		load_double(data, data_len, ans_col, row_idx);
	    break;
	}
	return;
}

static void load_tagval(const char *tag, int tag_len,
		const char *val, int val_len,
		SEXP ans, int row_idx)
{
	SEXP ans_names, col_name, ans_col;
	int ans_ncol0, j;

	ans_ncol0 = INTEGER(GET_ATTR(ans, install("ncol0")))[0];
	/* Lookup the current tag in 'names(ans)', starting from the end.
	   Since the number of tags in 'names(ans)' is typically very small
	   (< 25), we do this in a naive way (no hashing), because it's simple
	   and seems to be fast enough. */
	ans_names = GET_NAMES(ans);
	for (j = LENGTH(ans_names) - 1; j >= ans_ncol0; j--) {
		col_name = STRING_ELT(ans_names, j);
		if (LENGTH(col_name) == tag_len
		 && memcmp(CHAR(col_name), tag, tag_len) == 0)
			break;
	}
	if (j < ans_ncol0)
		return;  /* 'tag' was not found ==> nothing to do */
	ans_col = VECTOR_ELT(ans, j);
	load_string(val, val_len, ans_col, row_idx);
	return;
}

#define	IOBUF_SIZE 65536
static char errmsg_buf[200];

static void add_tag_to_buf(const char *tag, int tag_len, CharAEAE *tags_buf)
{
	CharAE **ae_p, *ae;
	int ntag, i;

	/* We want to store unique tags in 'tags_buf' so we first check
	   to see if 'tag' is already stored and don't do anything if it is.
	   Since the number of unique tags in a GFF file is typically very
	   small (< 25), we do this in a naive way (no hashing), because it's
	   simple and seems to be fast enough. */
	ntag = CharAEAE_get_nelt(tags_buf);
	for (i = 0, ae_p = tags_buf->elts; i < ntag; i++, ae_p++) {
		ae = *ae_p;
		if (CharAE_get_nelt(ae) == tag_len
		 && memcmp(ae->elts, tag, tag_len) == 0)
			return;  /* 'tag' was found ==> nothing to do */
	}
	/* 'tag' was not found ==> add it */
	ae = new_CharAE(tag_len);
	CharAE_set_nelt(ae, tag_len);
	memcpy(ae->elts, tag, tag_len);
	CharAEAE_insert_at(tags_buf, ntag, ae);
	return;
}

#define	UNKNOWN_FMT 0
#define	GFF2_FMT 2
#define	GFF3_FMT 3

/*
 * The tag-val components (i.e. the chunks between ';') with only white-space
 * characters are uninformative so we skip them. If all tag-val components have
 * only white-space characters then we return UNKNOWN_FMT. Otherwise the first
 * tag-val chunk with a non white-space character is used to determine the
 * format of the attributes.
 */
static int detect_attrcol_fmt(const char *data, int data_len)
{
	int only_spaces, i;
	char c;

	only_spaces = 1;
	for (i = 0; i < data_len; i++) {
		c = data[i];
		if (c == '=')
			return GFF3_FMT;
		if (c == '"')
			return GFF2_FMT;
		if (c == ';') {
			if (only_spaces)
				continue;
			return GFF2_FMT;
		}
		if (!isspace(c))
			only_spaces = 0;
	}
	return only_spaces ? UNKNOWN_FMT : GFF2_FMT;
}

static void parse_GFF3_tagval(const char *tagval, int tagval_len,
		SEXP ans, int row_idx, CharAEAE *tags_buf)
{
	int tag_len, val_len;
	char c;
	const char *val;

	for (tag_len = 0; tag_len < tagval_len; tag_len++) {
		c = tagval[tag_len];
		if (c == '=')
			break;
	}
	/* If 'tagval' is not in the tag=value format (i.e. if it has no =)
	   then we ignore it. */
	if (tag_len == tagval_len)
		return;
	if (ans != R_NilValue) {
		val = tagval + tag_len + 1;
		val_len = tagval_len - tag_len - 1;
		load_tagval(tagval, tag_len, val, val_len, ans, row_idx);
	}
	if (tags_buf != NULL)
		add_tag_to_buf(tagval, tag_len, tags_buf);
	return;
}

/*
 * It seems that embedded double-quotes might be allowed in the value part of
 * the tag value pairs of a GFF2 file, and that they are represented with 2
 * consecutive double-quotes. To handle them properly, we should replace the
 * 2 double-quotes by only 1, but this would require to generate a (shrinked)
 * copy of 'val'. So for now we don't do this and just leave the 2 consecutive
 * double-quotes in 'val' as-is, and attach an attribute to 'ans' to indicate
 * that we've found embedded double-quotes. We also issue a warning. If we ever
 * get that warning on real-world files, then we'll revisit this.
 */
static void check_for_embedded_dblquotes(const char *val, int val_len,
		SEXP ans)
{
	SEXP has_embedded_quotes;
	int i, j;

	has_embedded_quotes = GET_ATTR(ans, install("has_embedded_quotes"));
	if (!isNull(has_embedded_quotes) && LOGICAL(has_embedded_quotes)[0])
		return;
	for (i = 0, j = 1; j < val_len; i++, j++) {
		if (val[i] == '"' && val[j] == '"')
			break;
	}
	if (j == val_len)
		return;
	PROTECT(has_embedded_quotes = ScalarLogical(1));
	SET_ATTR(ans, install("has_embedded_quotes"), has_embedded_quotes);
	UNPROTECT(1);
	warning("the value part of some of the tag value pairs "
		"contains embedded double-quotes");
	return;
}

static void parse_GFF2_tagval(const char *tagval, int tagval_len,
		SEXP ans, int row_idx, CharAEAE *tags_buf)
{
	int i, tag_len, val_len;
	char c;
	const char *val;

	/* Trim leading space. */
	for (i = 0; i < tagval_len; i++) {
		c = tagval[i];
		if (!isspace(c))
			break;
	}
	tagval += i;
	tagval_len -= i;
	/* If 'tagval' has only white-space characters then we ignore it. */
	if (tagval_len == 0)
		return;
	for (tag_len = 0; tag_len < tagval_len; tag_len++) {
		c = tagval[tag_len];
		if (isspace(c))
			break;
	}
	if (ans != R_NilValue) {
		val = tagval + tag_len + 1;
		val_len = tagval_len - tag_len - 1;
		/* Trim leading space in 'val'. */
		for (i = 0; i < val_len; i++) {
			c = val[i];
			if (!isspace(c))
				break;
		}
		val += i;
		val_len -= i;
		/* Trim trailing space in 'val'. */
		for (i = val_len - 1; i >= 0; i--) {
			c = val[i];
			if (!isspace(c))
				break;
		}
		val_len = i + 1;
		/* Trim leading and trailing double-quotes in 'val'. */
		if (val_len != 0 && val[0] == '"') {
			val++;
			val_len--;
		}
		if (val_len != 0 && val[val_len - 1] == '"') {
			val_len--;
		}
		check_for_embedded_dblquotes(val, val_len, ans);
		load_tagval(tagval, tag_len, val, val_len, ans, row_idx);
	}
	if (tags_buf != NULL)
		add_tag_to_buf(tagval, tag_len, tags_buf);
	return;
}

static void parse_GFF3_attrcol(const char *data, int data_len,
		SEXP ans, int row_idx, CharAEAE *tags_buf)
{
	const char *tagval;
	int tagval_len, i;
	char c;

	tagval = data;
	tagval_len = 0;
	for (i = 0; i < data_len; i++) {
		c = data[i];
		if (c != ';') {
			tagval_len++;
			continue;
		}
		parse_GFF3_tagval(tagval, tagval_len, ans, row_idx, tags_buf);
		tagval = data + i + 1;
		tagval_len = 0;
	}
	parse_GFF3_tagval(tagval, tagval_len, ans, row_idx, tags_buf);
	return;
}

static void parse_GFF2_attrcol(const char *data, int data_len,
		SEXP ans, int row_idx, CharAEAE *tags_buf)
{
	const char *tagval;
	int tagval_len, in_quotes, i;
	char c;

	tagval = data;
	tagval_len = 0;
	in_quotes = 0;
	for (i = 0; i < data_len; i++) {
		c = data[i];
		if (c == '"') {
			tagval_len++;
			in_quotes = !in_quotes;
			continue;
		}
		if (in_quotes || c != ';') {
			tagval_len++;
			continue;
		}
		parse_GFF2_tagval(tagval, tagval_len, ans, row_idx, tags_buf);
		tagval = data + i + 1;
		tagval_len = 0;
	}
	parse_GFF2_tagval(tagval, tagval_len, ans, row_idx, tags_buf);
	return;
}

static const char *set_data_holders(Chars_holder *data_holders,
		const char *line, int lineno)
{
	int col_idx, i, data_len;
	const char *data;
	char c;

	col_idx = i = 0;
	data = line;
	data_len = 0;
	while ((c = line[i++])) {
		if (c != '\t') {
			data_len++;
			continue;
		}
		if (col_idx >= GFF_NCOL - 1)
			break;  /* 'c' contains the 9th tab */
		data_holders[col_idx].ptr = data;
		data_holders[col_idx].length = data_len;
		col_idx++;
		data = line + i;
		data_len = 0;
	}
	if (c) {
		/* We've seen 9 tabs but it's OK if the 9th tab is followed
		   by white spaces only (some GTF files are like that).
		   Otherwise we raise an error. */
		while ((c = line[i++])) {
			if (isspace(c))
				continue;
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "line %d has more than %d "
				 "tab-separated columns",
					 lineno, GFF_NCOL);
			return errmsg_buf;
		}
	}
	if (col_idx < GFF_NCOL - 2) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "line %d has less than %d tab-separated columns",
			 lineno, GFF_NCOL - 1);
		return errmsg_buf;
	}
	data_len = delete_trailing_LF_or_CRLF(data, data_len);
	data_holders[col_idx].ptr = data;
	data_holders[col_idx].length = data_len;
	return NULL;
}

static void check_filter(SEXP filter)
{
	int col_idx, nval, i;
	SEXP filter_elt, val;

	if (isNull(filter))
		return;
	if (!(IS_LIST(filter) && LENGTH(filter) == GFF_NCOL - 1))
		error("incorrect 'filter'");
	for (col_idx = 0; col_idx < GFF_NCOL - 1; col_idx++) {
		filter_elt = VECTOR_ELT(filter, col_idx);
		if (isNull(filter_elt))
			continue;
		if (!IS_CHARACTER(filter_elt))
			error("each list element in 'filter' must be "
			      "NULL or a character vector");
		nval = LENGTH(filter_elt);
		for (i = 0; i < nval; i++) {
			val = STRING_ELT(filter_elt, i);
			if (val == NA_STRING)
				error("'filter' cannot contain NAs");
		}
	}
	return;
}

static int pass_filter(Chars_holder *data_holders, SEXP filter)
{
	int col_idx, data_len, nval, i;
	SEXP filter_elt, val;
	const char *data;

	for (col_idx = 0; col_idx < GFF_NCOL - 1; col_idx++) {
		filter_elt = VECTOR_ELT(filter, col_idx);
		if (isNull(filter_elt))
			continue;
		data = data_holders[col_idx].ptr;
		data_len = data_holders[col_idx].length;
		nval = LENGTH(filter_elt);
		for (i = 0; i < nval; i++) {
			val = STRING_ELT(filter_elt, i);
			if (LENGTH(val) == data_len
			 && memcmp(CHAR(val), data, data_len) == 0)
				break;
		}
		if (i >= nval)
			return 0;
	}
	return 1;
}

static const char *parse_GFF_line(const char *line, int lineno,
		SEXP filter,
		SEXP ans, int *row_idx, const int *colmap0,
		int *attrcol_fmt,
		CharAEAE *tags_buf)
{
	Chars_holder data_holders[GFF_NCOL];
	const char *errmsg;
	int col_idx, data_len;
	const char *data;

	errmsg = set_data_holders(data_holders, line, lineno);
	if (errmsg != NULL)
		return errmsg;
	if (!(isNull(filter) || pass_filter(data_holders, filter)))
		return NULL;
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		data = data_holders[col_idx].ptr;
		data_len = data_holders[col_idx].length;
		if (ans != R_NilValue && colmap0[col_idx] != NA_INTEGER)
			load_data(data, data_len,
				  ans, *row_idx, col_idx, colmap0);
		if (col_idx != ATTRIBUTES_IDX)
			continue;
		if (*attrcol_fmt == UNKNOWN_FMT)
			*attrcol_fmt = detect_attrcol_fmt(data, data_len);
		if (*attrcol_fmt == GFF3_FMT)
			parse_GFF3_attrcol(data, data_len,
					   ans, *row_idx, tags_buf);
		else if (*attrcol_fmt == GFF2_FMT)
			parse_GFF2_attrcol(data, data_len,
					   ans, *row_idx, tags_buf);
	}
	(*row_idx)++;
	return NULL;
}

static const char *parse_GFF_file(SEXP filexp, SEXP filter,
		CharAEAE *tags_buf,	/* used during scan (1st pass) */
		int *ans_nrow,		/* used during scan (1st pass) */
		int *attrcol_fmt,	/* used during scan (1st pass) */
		CharAEAE *pragmas_buf,	/* used during scan (1st pass) */
		SEXP ans,		/* used during load (2nd pass) */
		const int *colmap0)	/* used during load (2nd pass) */
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
		if (c == '#') {
			if (buf[1] != '#')
				continue;  /* skip human-readable comment */
			if (pragmas_buf != NULL)
				append_string_to_CharAEAE(pragmas_buf, buf);
			continue;
		}
		if (c == '>')
			break;  /* stop parsing at first FASTA header */
		errmsg = parse_GFF_line(buf, lineno, filter,
					ans, &row_idx, colmap0,
					attrcol_fmt, tags_buf);
		if (errmsg != NULL)
			return errmsg;
	}
	*ans_nrow = row_idx;
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Performs the 1st pass of readGFF().
 * Args:
 *   filexp: A "file external pointer" (see src/io_utils.c in the XVector
 *           package).
 *   tags:   NULL or a character vector indicating which tags to load. NULL
 *           means load all tags.
 *   filter: NULL or a list of length GFF_NCOL-1. Each list element must be
 *           NULL or a character vector with no NAs.
 */
SEXP scan_gff(SEXP filexp, SEXP tags, SEXP filter)
{
	int ans_nrow, attrcol_fmt;
	CharAEAE *tags_buf, *pragmas_buf;
	const char *errmsg;
	SEXP scan_ans, scan_ans_elt;

	//init_clock("scan_gff: T1 = ");
	check_filter(filter);
	if (tags == R_NilValue) {
		tags_buf = new_CharAEAE(0, 0);
	} else {
		tags_buf = NULL;
	}
	attrcol_fmt = UNKNOWN_FMT;
	pragmas_buf = new_CharAEAE(0, 0);
	filexp_rewind(filexp);
	errmsg = parse_GFF_file(filexp, filter,
				tags_buf, &ans_nrow, &attrcol_fmt, pragmas_buf,
				R_NilValue, NULL);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);

	PROTECT(scan_ans = NEW_LIST(4));

	if (tags == R_NilValue) {
		PROTECT(scan_ans_elt = new_CHARACTER_from_CharAEAE(tags_buf));
		SET_ELEMENT(scan_ans, 0, scan_ans_elt);
		UNPROTECT(1);
	}
	PROTECT(scan_ans_elt = ScalarInteger(ans_nrow));
	SET_ELEMENT(scan_ans, 1, scan_ans_elt);
	UNPROTECT(1);
	PROTECT(scan_ans_elt = ScalarInteger(attrcol_fmt));
	SET_ELEMENT(scan_ans, 2, scan_ans_elt);
	UNPROTECT(1);
	PROTECT(scan_ans_elt = new_CHARACTER_from_CharAEAE(pragmas_buf));
	SET_ELEMENT(scan_ans, 3, scan_ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	//print_elapsed_time();
	return scan_ans;
}

/* --- .Call ENTRY POINT ---
 * Performs the 2nd pass of readGFF().
 * Args:
 *   filexp:      MUST be the same that was passed to scan_gff().
 *   tags:        A character vector indicating which tags to load. Unlike for
 *                scan_gff() above, it cannot be NULL.
 *   filter:      MUST be the same that was passed to scan_gff(). This time we
 *                don't check it.
 *   ans_nrow:    Integer vector of length 1 containing the number of rows
 *                scan_gff() said will end up in 'ans'.
 *   attrcol_fmt: 
 *   pragmas:     Character vector containing the pragma lines collected by
 *                scan_gff().
 *   colmap:      An integer vector of length GFF_NCOL indicating which columns
 *                to load and in which order.
 *   raw_data:    TRUE or FALSE. If TRUE, numeric columns (e.g. "start" or
 *                "score") are loaded as character vectors and as-is i.e. how
 *                they are found in the file.
 */
SEXP load_gff(SEXP filexp, SEXP tags, SEXP filter,
		SEXP ans_nrow, SEXP attrcol_fmt, SEXP pragmas,
		SEXP colmap, SEXP raw_data)
{
	int colmap0[GFF_NCOL], ans_ncol0;
	SEXP ans;
	const char *errmsg;

	//init_clock("load_gff: T2 = ");
	ans_ncol0 = prepare_colmap0(colmap0, colmap);
	PROTECT(ans = alloc_ans(INTEGER(ans_nrow)[0], ans_ncol0, colmap0,
				tags, pragmas, raw_data));
	filexp_rewind(filexp);
	errmsg = parse_GFF_file(filexp, filter,
				NULL, INTEGER(ans_nrow),
				INTEGER(attrcol_fmt), NULL,
				ans, colmap0);
	UNPROTECT(1);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);
	//print_elapsed_time();
	return ans;
}

