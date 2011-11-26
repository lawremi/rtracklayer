#include "utils.h"

SEXP _STRSXP_collapse(SEXP x, SEXP sep) {
  SEXP ans;
  int len = 0;
  char *collapsed, *dest;
  char c_sep = CHAR(STRING_ELT(sep, 0))[0];
  
  if (TYPEOF(x) != STRSXP)
    error("_STRSXP_collapse: expected a STRSXP");
  if (length(x) == 1)
    return(STRING_ELT(x, 0));
  
  for (int i = 0; i < length(x); i++)
    len += strlen(CHAR(STRING_ELT(x, i))) + 1;
  collapsed = R_alloc(sizeof(char), len);

  dest = collapsed;
  for (int i = 0; i < length(x); i++) {
    const char *src = CHAR(STRING_ELT(x, i));
    int src_len = strlen(src);
    strcpy(dest, src);
    dest[src_len] = c_sep;
    dest += src_len + 1;
  }

  return mkCharLen(collapsed, len - 1 * (length(x) > 0));
}

/* Utility for collapsing elements of a CharacterList */
/* Assumes that 'x' is a 'list' */
SEXP CharacterList_pasteCollapse(SEXP x, SEXP sep) {
  SEXP ans;
  if (TYPEOF(x) != VECSXP)
    error("CharacterList_collapse: expected a list");
  PROTECT(ans = allocVector(STRSXP, length(x)));
  for (int i = 0; i < length(x); i++)
    SET_STRING_ELT(ans, i, _STRSXP_collapse(VECTOR_ELT(x, i), sep));
  UNPROTECT(1);
  return ans;
}
