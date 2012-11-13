### =========================================================================
### Utilities
### -------------------------------------------------------------------------

pasteCollapse <- function(x, collapse = ",") {
  if (!is(x, "CharacterList"))
    stop("'x' must be a CharacterList")
  if (!isSingleString(collapse) || nchar(collapse) != 1L)
    stop("'collapse' must be a single string, with a single character")
  x <- as.list(x)
  .Call(CharacterList_pasteCollapse, x, collapse)
}

asRangedData.warning.msg <- function(fname, if.FALSE="GRanges",
                                            if.TRUE="RangedData") {
  msg <- c("Starting with BioC 2.12, the default value for the ",
           "'asRangedData' argument\n  of %s() has changed ",
           "from TRUE to FALSE. ",
           "Explicitly provide a\n  value for this argument to avoid ",
           "this warning e.g. with\n",
           "    %s(..., asRangedData=FALSE)\n",
           "  to get a %s object, or with\n",
           "    %s(..., asRangedData=TRUE)\n",
           "  to get a %s object. ",
           "This warning will be removed in BioC 2.13")
  fmt <- paste0(msg, collapse="")
  sprintf(fmt, fname, fname, if.FALSE, fname, if.TRUE)
}

