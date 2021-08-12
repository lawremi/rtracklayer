#include "XVector_interface.h"
#include "S4Vectors_interface.h"

#include <R_ext/Connections.h>  /* for R_ReadConnection() */

#include <ctype.h>   /* for isspace() */
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


/****************************************************************************
 * filexp_gets2(): A version of filexp_gets() that also works on connections
 */

Rconnection getConnection(int n);  /* not in <R_ext/Connections.h>, why? */

static char con_buf[25000];
static int con_buf_len, con_buf_offset;

static void init_con_buf()
{
	con_buf_len = con_buf_offset = 0;
	return;
}

static int filexp_gets2(SEXP filexp, char *buf, int buf_size, int *EOL_in_buf)
{
	Rconnection con;
	int buf_offset;
	char c;

	if (TYPEOF(filexp) == EXTPTRSXP)
		return filexp_gets(filexp, buf, buf_size, EOL_in_buf);
	buf_offset = *EOL_in_buf = 0;
	while (buf_offset < buf_size - 1) {
		if (con_buf_offset == con_buf_len) {
			con = getConnection(asInteger(filexp));
			con_buf_len = (int) R_ReadConnection(con,
					con_buf,
					sizeof(con_buf) / sizeof(char));
			if (con_buf_len == 0)
				break;
			con_buf_offset = 0;
		}
		c = con_buf[con_buf_offset++];
		buf[buf_offset++] = c;
		if (c == '\n') {
			*EOL_in_buf = 1;
			break;
		}
	}
	buf[buf_offset] = '\0';
	if (buf_offset == 0)
		return 0;
	if (con_buf_len == 0 || *EOL_in_buf)
		return 2;
	return 1;
}


/****************************************************************************
 * Buffer of tags
 */

#define	MAX_COLLECTED_TAGS 4096
typedef struct tags_buf {
	CharAEAE *collected_tags;  /* used in COLLECT MODE */
	SEXP tags;		   /* used in LOAD MODE */
	/* We use a hash table for fast mapping of the tags to their
	   corresponding column in the output of readGFF(). */
	struct htab htab;
} TagsBuf;

#define	IN_COLLECT_MODE(tags_buf) ((tags_buf)->collected_tags != NULL)

/* Copied from R_HOME/src/main/envir.c */
static unsigned int char_hash(const char *s, int len)
{
	/* djb2 as from http://www.cse.yorku.ca/~oz/hash.html */
	int i;
	unsigned int h = 5381;

	for (i = 0; i < len; i++)
		h = ((h << 5) + h) + (*s);
	return h;
}

static int TagsBuf_get_bucket_idx(const TagsBuf *tags_buf,
		const char *tag, int tag_len)
{
	unsigned int hval;
	int bucket_idx, M, i2;
	const int *buckets;
	const CharAE *tag_ae;
	SEXP tags_elt;

	hval = char_hash(tag, tag_len);
	bucket_idx = hval & tags_buf->htab.Mminus1;
	buckets = tags_buf->htab.buckets;
	M = tags_buf->htab.M;
	while ((i2 = buckets[bucket_idx]) != NA_INTEGER) {
		if (IN_COLLECT_MODE(tags_buf)) {
			/* COLLECT MODE */
			tag_ae = tags_buf->collected_tags->elts[i2];
			if (CharAE_get_nelt(tag_ae) == tag_len
			 && memcmp(tag_ae->elts, tag, tag_len) == 0)
				break;
		} else {
			/* LOAD MODE */
			tags_elt = STRING_ELT(tags_buf->tags, i2);
			if (LENGTH(tags_elt) == tag_len
			 && memcmp(CHAR(tags_elt), tag, tag_len) == 0)
				break;
		}
		bucket_idx = (bucket_idx + 1) % M;
	}
	return bucket_idx;
}

static int TagsBuf_get_bucket_val(const TagsBuf *tags_buf, int bucket_idx)
{
	return get_hbucket_val(&(tags_buf->htab), bucket_idx);
}

static void TagsBuf_set_bucket_val(TagsBuf *tags_buf, int bucket_idx, int val)
{
	set_hbucket_val(&(tags_buf->htab), bucket_idx, val);
	return;
}

static void init_tags_buf_in_COLLECT_MODE(TagsBuf *tags_buf)
{
	tags_buf->collected_tags = new_CharAEAE(MAX_COLLECTED_TAGS, 0);
	tags_buf->tags = R_NilValue;
	tags_buf->htab = new_htab(MAX_COLLECTED_TAGS);
	return;
}

/* The tags in 'tags' must be unique. */
static void init_tags_buf_in_LOAD_MODE(TagsBuf *tags_buf, SEXP tags)
{
	int ntag, i, tag_len, bucket_idx, bucket_val;
	SEXP tags_elt;
	const char *tag;

	tags_buf->collected_tags = NULL;
	tags_buf->tags = tags;
	ntag = LENGTH(tags);
	tags_buf->htab = new_htab(ntag);
	for (i = 0; i < ntag; i++) {
		tags_elt = STRING_ELT(tags, i);
		if (tags_elt == NA_STRING)
			error("'tags' cannot contain NAs");
		tag = CHAR(tags_elt);
		tag_len = LENGTH(tags_elt);
		bucket_idx = TagsBuf_get_bucket_idx(tags_buf, tag, tag_len);
		bucket_val = TagsBuf_get_bucket_val(tags_buf, bucket_idx);
		if (bucket_val != NA_INTEGER)
			error("'tags' cannot contain duplicates");
		TagsBuf_set_bucket_val(tags_buf, bucket_idx, i);
	}
	return;
}

static void collect_tag(TagsBuf *tags_buf, const char *tag, int tag_len)
{
	int bucket_idx, bucket_val;
	int i;
	CharAE *ae;

	/* We want to store unique tags in 'tags_buf' so we first check to see
	   if 'tag' is already stored and we don't do anything if it is. */
	bucket_idx = TagsBuf_get_bucket_idx(tags_buf, tag, tag_len);
	bucket_val = TagsBuf_get_bucket_val(tags_buf, bucket_idx);
	if (bucket_val != NA_INTEGER)
		return;  /* 'tag' was found ==> nothing to do */
	/* 'tag' was not found ==> add it */
	i = CharAEAE_get_nelt(tags_buf->collected_tags);
	if (i >= MAX_COLLECTED_TAGS)
		error("GFF files with more than %d tags are not supported",
		      MAX_COLLECTED_TAGS);
	TagsBuf_set_bucket_val(tags_buf, bucket_idx, i);
	ae = new_CharAE(tag_len);
	CharAE_set_nelt(ae, tag_len);
	memcpy(ae->elts, tag, tag_len);
	CharAEAE_insert_at(tags_buf->collected_tags, i, ae);
	return;
}

static int match_tag(const TagsBuf *tags_buf, const char *tag, int tag_len)
{
	int bucket_idx;

	bucket_idx = TagsBuf_get_bucket_idx(tags_buf, tag, tag_len);
	return TagsBuf_get_bucket_val(tags_buf, bucket_idx);
}

static SEXP get_collected_tags(const TagsBuf *tags_buf)
{
	if (tags_buf == NULL || !IN_COLLECT_MODE(tags_buf))
		return R_NilValue;
	return new_CHARACTER_from_CharAEAE(tags_buf->collected_tags);
}


/****************************************************************************
 * gff_colnames()
 *
 * See http://www.sequenceontology.org/resources/gff3.html for the official
 * GFF3 specs.
 */

static const char *col_names[] = {
	"seqid",
	"source",
	"type",
	"start",
	"end",
	"score",
	"strand",
	"phase",
	"attributes"  /* "group" for GFF1 */
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

#define	GFF_NCOL ((int) (sizeof(col_names) / sizeof(char *)))

#define	SEQID_IDX 0
#define	SOURCE_IDX 1
#define	TYPE_IDX 2
#define	START_IDX 3
#define	END_IDX 4
#define	SCORE_IDX 5
#define	STRAND_IDX 6
#define	PHASE_IDX 7
#define	ATTRIBUTES_IDX 8

static const char *gff_colname(int col_idx, int gff1)
{
	if (col_idx == ATTRIBUTES_IDX && gff1)
		return "group";
	return col_names[col_idx];
}

/* --- .Call ENTRY POINT --- */
SEXP gff_colnames(SEXP GFF1)
{
	SEXP ans, ans_elt;
	int col_idx;
	const char *colname;

	PROTECT(ans = NEW_CHARACTER(GFF_NCOL));
	for (col_idx = 0; col_idx < GFF_NCOL; col_idx++) {
		colname = gff_colname(col_idx, LOGICAL(GFF1)[0]);
		PROTECT(ans_elt = mkChar(colname));
		SET_STRING_ELT(ans, col_idx, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * read_gff_pragmas(), scan_gff(), and load_gff()
 */

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
		SEXP tags, SEXP pragmas, SEXP attrcol_fmt, SEXP raw_data)
{
	int ans_ntag, ans_ncol, gff1, is_raw, col_idx, j, i;
	SEXP ans, ans_attr, ans_names, ans_col, ans_colname, tags_elt;
	SEXPTYPE col_type;
	const char *colname;

	ans_ntag = LENGTH(tags);
	ans_ncol = ans_ncol0 + ans_ntag;
	gff1 = INTEGER(attrcol_fmt)[0] == 1;
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
		colname = gff_colname(col_idx, gff1);
		PROTECT(ans_colname = mkChar(colname));
		SET_STRING_ELT(ans_names, j, ans_colname);
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
		PROTECT(ans_colname = duplicate(tags_elt));
		SET_STRING_ELT(ans_names, j, ans_colname);
		UNPROTECT(1);
	}

	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* list_as_data_frame() performs IN-PLACE coercion */
	list_as_data_frame(ans, ans_nrow);

	/* Set additional attributes. */
	PROTECT(ans_attr = duplicate(pragmas));
	SET_ATTR(ans, install("pragmas"), ans_attr);
	UNPROTECT(1);
	PROTECT(ans_attr = duplicate(attrcol_fmt));
	SET_ATTR(ans, install("attrcol_fmt"), ans_attr);
	UNPROTECT(1);
	PROTECT(ans_attr = ScalarInteger(ans_ncol0));
	SET_ATTR(ans, install("ncol0"), ans_attr);
	UNPROTECT(1);
	PROTECT(ans_attr = ScalarInteger(ans_ntag));
	SET_ATTR(ans, install("ntag"), ans_attr);
	UNPROTECT(1);
	PROTECT(ans_attr = duplicate(raw_data));
	SET_ATTR(ans, install("raw_data"), raw_data);
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

static char errmsg_buf[256];

static const char *load_data(const char *data, int data_len,
		SEXP ans, int row_idx, int col_idx, const int *colmap0,
		int lineno)
{
	SEXP ans_col;
	int is_raw;
	SEXPTYPE col_type;

	ans_col = VECTOR_ELT(ans, colmap0[col_idx]);
	is_raw = LOGICAL(GET_ATTR(ans, install("raw_data")))[0];
	if (is_raw) {
		load_string(data, data_len, ans_col, row_idx);
		return NULL;
	}
	col_type = col_types[col_idx];
	if (col_idx == START_IDX || col_idx == END_IDX) {
		reset_ovflow_flag();
		load_int(data, data_len, ans_col, row_idx);
		if (get_ovflow_flag()) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "line %d contains values greater than 2^31-1 "
				 "\n  (= .Machine$integer.max) in column 4 "
				 "(start) and/or 5 (end).\n  Bioconductor "
				 "does not support such GFF files at the "
				 "moment. Sorry!", lineno);
			return errmsg_buf;
		}
		return NULL;
	}
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
	return NULL;
}

static void load_tagval(const char *tag, int tag_len,
		const char *val, int val_len,
		SEXP ans, int row_idx, const TagsBuf *tags_buf)
{
	int j;
	SEXP ans_col;

	j = match_tag(tags_buf, tag, tag_len);
	if (j == NA_INTEGER)
		return;  /* 'tag' was not found ==> nothing to do */
	j += INTEGER(GET_ATTR(ans, install("ncol0")))[0];
	ans_col = VECTOR_ELT(ans, j);
	load_string(val, val_len, ans_col, row_idx);
	return;
}

/* Some stupid GFF files like GCF_000001735.4_TAIR10.1_genomic.gff.gz (RefSeq
 * GFF for Arabidopsis thaliana TAIR10.1) can contain lines up to 168975
 * chars long. */
#define	IOBUF_SIZE 200000

/*
 * We use a heuristic to detect the format of the "attributes" col.
 * Terminology:
 *   - Chunks in 'data' separated by ';' are called units.
 *   - A unit with only white-space characters is called a "white" unit.
 *   - A "tag-like word" is a sequence of contiguous non-white-space characters
 *     with no '=' in it.
 *   - A unit made of one "tag-like word" (possibly surrounded by white-space
 *     characters) is called a "tag-only" unit.
 * Heuristic:
 *   (a) If 'data' contains only 1 unit (i.e. no ';' in it),
 *       then:
 *         - if it's "white" -> return UNKNOWN_FMT
 *         - if it's "tag-only" -> return GFF1_FMT
 *         - otherwise -> (c) below applies.
 *   (b) If 'data' contains > 1 units (i.e. at least 1 ';' in it), then the
 *       format can't be GFF1_FMT anymore so we have to choose between GFF2_FMT
 *       and GFF3_FMT. For that matter "white" and "tag-only" units are
 *       considered uninformative so we skip them. If all units are
 *       uninformative then we return UNKNOWN_FMT. Otherwise the first
 *       informative unit is used and (c) below applies.
 *   (c) Now we're looking at a unit that is not "white" or "tag-only" and
 *       want to be able to tell whether its format is GFF2_FMT or GFF3_FMT
 *       (UNKNOWN_FMT or GFF1_FMT are not an option anymore). The rule is: if
 *       the unit contains a '=' and if this '=' is preceded by one
 *       "tag-like word" only (possibly surrounded by white-space characters)
 *       then the format is GFF3_FMT. Otherwise it's GFF2_FMT.
 * Examples:
 *    data                   detected format
 *   |ID=4|----------------> GFF3_FMT
 *   |  ID  = 4   |--------> GFF3_FMT
 *   |=|-------------------> GFF3_FMT
 *   |     = 4   |---------> GFF3_FMT
 *   |   = ID 4    |-------> GFF3_FMT
 *   |ID 4|----------------> GFF2_FMT
 *   |  ID   4  |----------> GFF2_FMT
 *   |X Y=4|---------------> inherently ambiguous (could be GFF2 or GFF3) but
 *                           our heuristic returns GFF2_FMT (tag is "X" and
 *                           value is "Y=4")
 *   |XY z =|--------------> GFF2_FMT (tag is "XY" and value is "z =")
 *   |; ;ID=4;Parent 12|---> GFF3_FMT (we look at the 1st informative unit)
 *   |; ;ID 4;Parent=12|---> GFF2_FMT (we look at the 1st informative unit)
 *   |99|------------------> GFF1_FMT
 *   ||--------------------> UNKNOWN_FMT
 *   |   |-----------------> UNKNOWN_FMT
 *   | ; 99 ;   |----------> UNKNOWN_FMT
 *   |99; ID =  4|---------> GFF3_FMT
 */
#define	UNKNOWN_FMT 0
#define	GFF1_FMT 1
#define	GFF2_FMT 2
#define	GFF3_FMT 3

#define	FIRST_SPACE 0
#define	FIRST_WORD 1
#define	SECOND_SPACE 2
static int detect_attrcol_fmt(const char *data, int data_len)
{
	int nsep, state, i;
	char c;

	nsep = 0;
	state = FIRST_SPACE;
	for (i = 0; i < data_len; i++) {
		c = data[i];
		if (isspace(c)) {
			if (state == FIRST_WORD)
				state = SECOND_SPACE;
			continue;
		}
		if (c == ';') {
			nsep++;
			/* We came to the end of a unit that was "white" (i.e.
			   current state is FIRST_SPACE) or "tag-only" (i.e.
			   current state is FIRST_WORD or SECOND_SPACE). This
			   is considered uninformative i.e. it didn't allow us
			   to decide between GFF2_FMT and GFF3_FMT (see (b)
			   above). */
			state = FIRST_SPACE;
			continue;
		}
		if (c == '=')
			return GFF3_FMT;
		if (state == SECOND_SPACE)
			return GFF2_FMT;
		if (state == FIRST_SPACE)
			state = FIRST_WORD;
	}
	if (nsep == 0 && state != FIRST_SPACE)
		return GFF1_FMT;
	return UNKNOWN_FMT;
}

static void parse_GFF3_tagval(const char *tagval, int tagval_len,
		SEXP ans, int row_idx, TagsBuf *tags_buf)
{
	int tag_len, val_len;
	char c;
	const char *val;

	/* Some GFF3 files have a space betwwen the tag=value and the
	   preceding ; e.g.
	       ID=Csa1M000010.1; Parent=Csa1G000010; Name=Csa1M000010.1
	   We skip those spaces. */
	while (tagval_len > 0 && tagval[0] == ' ') {
		tagval++;
		tagval_len--;
	}
	/* Compute 'tag_len'. */
	for (tag_len = 0; tag_len < tagval_len; tag_len++) {
		c = tagval[tag_len];
		if (c == '=')
			break;
	}
	/* If 'tagval' is not in the "tag=value" format then we ignore it. */
	if (tag_len >= tagval_len)
		return;
	if (ans != R_NilValue) {
		val = tagval + tag_len + 1;
		val_len = tagval_len - tag_len - 1;
		load_tagval(tagval, tag_len, val, val_len,
			    ans, row_idx, tags_buf);
		return;
	}
	if (tags_buf != NULL && IN_COLLECT_MODE(tags_buf))
		collect_tag(tags_buf, tagval, tag_len);
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
	if (j >= val_len)
		return;
	PROTECT(has_embedded_quotes = ScalarLogical(1));
	SET_ATTR(ans, install("has_embedded_quotes"), has_embedded_quotes);
	UNPROTECT(1);
	warning("the value part of some of the tag value pairs "
		"contains embedded double-quotes");
	return;
}

static void parse_GFF2_tagval(const char *tagval, int tagval_len,
		SEXP ans, int row_idx, TagsBuf *tags_buf)
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
	/* Compute 'tag_len'. */
	for (tag_len = 0; tag_len < tagval_len; tag_len++) {
		c = tagval[tag_len];
		if (isspace(c))
			break;
	}
	/* If 'tagval' is not in the "tag value" format then we ignore it. */
	if (tag_len >= tagval_len)
		return;
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
		load_tagval(tagval, tag_len, val, val_len,
			    ans, row_idx, tags_buf);
		return;
	}
	if (tags_buf != NULL && IN_COLLECT_MODE(tags_buf))
		collect_tag(tags_buf, tagval, tag_len);
	return;
}

static void parse_GFF3_attrcol(const char *data, int data_len,
		SEXP ans, int row_idx, TagsBuf *tags_buf)
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
		SEXP ans, int row_idx, TagsBuf *tags_buf)
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
		data_holders[col_idx].ptr = data;
		data_holders[col_idx].length = data_len;
		col_idx++;
		if (col_idx == GFF_NCOL)
			break;
		data = line + i;
		data_len = 0;
	}
	if (col_idx == GFF_NCOL) {
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
	} else {
		if (col_idx < GFF_NCOL - 2) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "line %d has less than %d "
				 "tab-separated columns",
				 lineno, GFF_NCOL - 1);
			return errmsg_buf;
		}
		data_len = delete_trailing_LF_or_CRLF(data, data_len);
		data_holders[col_idx].ptr = data;
		data_holders[col_idx].length = data_len;
		col_idx++;
		if (col_idx == GFF_NCOL - 1)
			data_holders[col_idx].length = 0;
	}
	return NULL;
}

static void check_filter(SEXP filter, int attrcol_fmt)
{
	int filter_len, col_idx, nval, i;
	SEXP filter_elt, val;

	if (isNull(filter))
		return;
	filter_len = attrcol_fmt == 1 ? GFF_NCOL : GFF_NCOL - 1;
	if (!(IS_LIST(filter) && LENGTH(filter) == filter_len))
		error("incorrect 'filter'");
	for (col_idx = 0; col_idx < filter_len; col_idx++) {
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
	int filter_len, col_idx, data_len, nval, i;
	SEXP filter_elt, val;
	const char *data;

	filter_len = LENGTH(filter);
	for (col_idx = 0; col_idx < filter_len; col_idx++) {
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
		int *attrcol_fmt,
		SEXP filter,		/* used during scan and load */
		SEXP ans,		/* used during load (2nd pass) */
		int *row_idx,		/* used during scan and load */
		const int *colmap0,	/* used during load (2nd pass) */
		TagsBuf *tags_buf)	/* used during scan and load */
{
	const char *errmsg;
	Chars_holder data_holders[GFF_NCOL];
	const Chars_holder *data_holder;
	const char *data;
	int data_len, col_idx;

	errmsg = set_data_holders(data_holders, line, lineno);
	if (errmsg != NULL)
		return errmsg;
	/* Try to detect the format of the "attributes" col before
	   filtering. */
	if (*attrcol_fmt == UNKNOWN_FMT) {
		data_holder = data_holders + ATTRIBUTES_IDX;
		data = data_holder->ptr;
		data_len = data_holder->length;
		*attrcol_fmt = detect_attrcol_fmt(data, data_len);
	}
	if (!(isNull(filter) || pass_filter(data_holders, filter)))
		return NULL;
	for (col_idx = 0, data_holder = data_holders;
	     col_idx < GFF_NCOL;
	     col_idx++, data_holder++)
	{
		data = data_holder->ptr;
		data_len = data_holder->length;
		if (ans != R_NilValue && colmap0[col_idx] != NA_INTEGER) {
			errmsg = load_data(data, data_len,
					   ans, *row_idx, col_idx, colmap0,
					   lineno);
			if (errmsg != NULL)
				return errmsg;
		}
		if (col_idx != ATTRIBUTES_IDX)
			continue;
		switch (*attrcol_fmt) {
		    case GFF3_FMT:
			parse_GFF3_attrcol(data, data_len,
					   ans, *row_idx, tags_buf);
		    break;
		    case GFF2_FMT:
			parse_GFF2_attrcol(data, data_len,
					   ans, *row_idx, tags_buf);
		    break;
		}
	}
	(*row_idx)++;
	return NULL;
}

/* Stop parsing after 1st GFF line. */
static const char *load_GFF_pragmas(SEXP filexp,
		CharAEAE *pragmas_buf, int *attrcol_fmt)
{
	int row_idx, lineno, ret_code, EOL_in_buf, buf_len;
	char buf[IOBUF_SIZE], c;

	if (TYPEOF(filexp) != EXTPTRSXP)
		init_con_buf();
	row_idx = 0;
	for (lineno = 1;
	     (ret_code = filexp_gets2(filexp, buf, IOBUF_SIZE, &EOL_in_buf));
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
		if (c != '#') {
			/* Skip empty line. */
			if (c == '\n' || (c == '\r' && buf[1] == '\n'))
				continue;
			/* Stop parsing at 1st FASTA header... */
			if (c == '>')
				return NULL;
			/* ... or at 1st GFF line. */
			return parse_GFF_line(buf, lineno, attrcol_fmt,
					      R_NilValue,
					      R_NilValue, &row_idx, NULL, NULL);
		}
		/* Line starting with a single # -> human-readable comment. */
		if (buf[1] != '#')
			continue;
		/* Line starting with ## -> pragma line. */
		buf_len = delete_trailing_LF_or_CRLF(buf, -1);
		buf[buf_len] = '\0';
		CharAEAE_append_string(pragmas_buf, buf);
	}
	return NULL;
}

static const char *parse_GFF_file(SEXP filexp, int *attrcol_fmt,
		SEXP filter,		/* used during scan and load */
		int *nrows,		/* used during scan and load */
		SEXP ans,		/* used during load (2nd pass) */
		const int *colmap0,	/* used during load (2nd pass) */
		TagsBuf *tags_buf)	/* used during scan and load */
{
	int row_idx, lineno, ret_code, EOL_in_buf;
	char buf[IOBUF_SIZE], c;
	const char *errmsg;

	if (TYPEOF(filexp) != EXTPTRSXP)
		init_con_buf();
	row_idx = 0;
	for (lineno = 1;
	     (ret_code = filexp_gets2(filexp, buf, IOBUF_SIZE, &EOL_in_buf));
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
		if (*nrows >= 0 && row_idx == *nrows)
			return NULL;
		c = buf[0];
		if (c == '#')
			continue;
		if (c == '\n' || (c == '\r' && buf[1] == '\n'))
			continue;  /* skip empty line */
		if (c == '>')
			break;  /* stop parsing at 1st FASTA header */
		errmsg = parse_GFF_line(buf, lineno, attrcol_fmt,
					filter,
					ans, &row_idx, colmap0, tags_buf);
		if (errmsg != NULL)
			return errmsg;
	}
	*nrows = row_idx;
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Return pragmas and detected format of the "attributes" col from the 1st
 * GFF line.
 * Arg:
 *   filexp: A "file external pointer" (see src/io_utils.c in the XVector
 *           package). Alternatively can be a connection.
 */
SEXP read_gff_pragmas(SEXP filexp)
{
	CharAEAE *pragmas_buf;
	int attrcol_fmt0;
	const char *errmsg;
	SEXP pragmas, attrcol_fmt;

	pragmas_buf = new_CharAEAE(0, 0);
	attrcol_fmt0 = UNKNOWN_FMT;
	errmsg = load_GFF_pragmas(filexp, pragmas_buf, &attrcol_fmt0);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);
	PROTECT(pragmas = new_CHARACTER_from_CharAEAE(pragmas_buf));
	PROTECT(attrcol_fmt = ScalarInteger(attrcol_fmt0));
	SET_ATTR(pragmas, install("attrcol_fmt"), attrcol_fmt);
	UNPROTECT(2);
	return pragmas;
}

/* --- .Call ENTRY POINT ---
 * Perform the 1st pass of readGFF().
 * Args:
 *   filexp:      MUST be the same as passed to read_gff_pragmas().
 *   attrcol_fmt:
 *   tags:        NULL or a character vector indicating which tags to load.
 *                NULL means load all tags.
 *   filter:      NULL or a list of length GFF_NCOL-1 (or GFF_NCOL if
 *                'attrcol_fmt' is 1). Each list element must be NULL or a
 *                character vector with no NAs.
 *   nrows:       -1 or the maximum number of rows to read in (after
 *                filtering).
 */
SEXP scan_gff(SEXP filexp, SEXP attrcol_fmt, SEXP tags,
	      SEXP filter, SEXP nrows)
{
	int attrcol_fmt0, nrows0;
	TagsBuf tags_buf, *tags_buf_p;
	const char *errmsg;
	SEXP scan_ans, scan_ans_elt;

	//init_clock("scan_gff: T1 = ");
	attrcol_fmt0 = INTEGER(attrcol_fmt)[0];
	if (tags == R_NilValue) {
		tags_buf_p = &tags_buf;
		init_tags_buf_in_COLLECT_MODE(tags_buf_p);
	} else {
		tags_buf_p = NULL;
	}
	check_filter(filter, attrcol_fmt0);
	nrows0 = INTEGER(nrows)[0];
	errmsg = parse_GFF_file(filexp, &attrcol_fmt0,
				filter, &nrows0,
				R_NilValue, NULL, tags_buf_p);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);

	PROTECT(scan_ans = NEW_LIST(2));
	PROTECT(scan_ans_elt = get_collected_tags(tags_buf_p));
	SET_ELEMENT(scan_ans, 0, scan_ans_elt);
	UNPROTECT(1);
	PROTECT(scan_ans_elt = ScalarInteger(nrows0));
	SET_ELEMENT(scan_ans, 1, scan_ans_elt);
	UNPROTECT(2);
	//print_elapsed_time();
	return scan_ans;
}

/* --- .Call ENTRY POINT ---
 * Perform the 2nd pass of readGFF().
 * Args:
 *   filexp:      MUST be the same as passed to scan_gff().
 *   attrcol_fmt:
 *   tags:        A character vector indicating which tags to load. Unlike for
 *                scan_gff() above, it cannot be NULL.
 *   filter:      MUST be the same as passed to scan_gff(). This time we
 *                don't check it.
 *   nrows:       Integer vector of length 1 containing the number of rows
 *                scan_gff() said will end up in 'ans'.
 *   pragmas:     Character vector containing the pragma lines returned by
 *                read_gff_pragmas().
 *   colmap:      An integer vector of length GFF_NCOL indicating which columns
 *                to load and in which order.
 *   raw_data:    TRUE or FALSE. If TRUE, numeric columns (e.g. "start" or
 *                "score") are loaded as character vectors and as-is i.e. how
 *                they are found in the file.
 */
SEXP load_gff(SEXP filexp, SEXP attrcol_fmt, SEXP tags,
	      SEXP filter, SEXP nrows,
	      SEXP pragmas, SEXP colmap, SEXP raw_data)
{
	int attrcol_fmt0, colmap0[GFF_NCOL], ans_ncol0;
	TagsBuf tags_buf;
	SEXP ans;
	const char *errmsg;

	//init_clock("load_gff: T2 = ");
	attrcol_fmt0 = INTEGER(attrcol_fmt)[0];
	init_tags_buf_in_LOAD_MODE(&tags_buf, tags);
	ans_ncol0 = prepare_colmap0(colmap0, colmap);
	PROTECT(ans = alloc_ans(INTEGER(nrows)[0], ans_ncol0, colmap0,
				tags, pragmas, attrcol_fmt, raw_data));
	errmsg = parse_GFF_file(filexp, &attrcol_fmt0,
				filter, INTEGER(nrows),
				ans, colmap0, &tags_buf);
	UNPROTECT(1);
	if (errmsg != NULL)
		error("reading GFF file: %s", errmsg);
	//print_elapsed_time();
	return ans;
}

