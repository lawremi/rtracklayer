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

.asRangedData_error_msg <- function(fname, if.FALSE="GRanges",
                                             if.TRUE="RangedData") {
  msg <- c("  Starting with BioC 2.14, calling %s() with 'asRangedData=TRUE' ",
           "is\n  unsupported. If you wish %s() to return a %s object,",
           "\n  then you can coerce the returned object with '",
           "as(..., \"%s\")'.\n",
           "  However, we strongly recommend that you start migrating your ",
           "code to\n  operate on %s objects instead of %s objects.\n",
           "  Please ask on the bioc-devel mailing list if you ",
           "have questions\n  or concerns about this ",
           "(http://bioconductor.org/help/mailing-list/)")
  fmt <- paste0(msg, collapse="")
  sprintf(fmt, fname, fname, if.TRUE, if.TRUE, if.FALSE, if.TRUE, if.TRUE)
}

normarg_asRangedData <- function(asRangedData, fname, if.FALSE="GRanges",
                                                      if.TRUE="RangedData") {
  if (!isTRUEorFALSE(asRangedData))
    stop("'asRangedData' must be TRUE or FALSE")
  if (asRangedData)
    stop(.asRangedData_error_msg(fname, if.FALSE, if.TRUE))
  asRangedData
}

