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

.asRangedData_warning_msg <- function(fname, if.FALSE="GRanges",
                                             if.TRUE="RangedData") {
  msg <- c("  Starting with BioC 2.13, calling %s() with 'asRangedData=TRUE' ",
           "is\n  deprecated. If you wish %s() to return a %s object,",
           "\n  then you can coerce the returned object with '",
           "as(..., \"%s\")'.\n",
           "  However, we strongly recommend that you start migrating your ",
           "code to\n  operate on %s objects instead of %s objects.\n",
           "  Starting with BioC 2.14, support for %s objects will be ",
           "limited.\n  Please ask on the bioc-devel mailing list if you ",
           "have questions\n  or concerns about this ",
           "(http://bioconductor.org/help/mailing-list/)")
  fmt <- paste0(msg, collapse="")
  sprintf(fmt, fname, fname, if.TRUE, if.TRUE, if.FALSE, if.TRUE, if.TRUE)
}

normarg_asRangedData <- function(asRangedData, fname, if.FALSE="GRanges",
                                                      if.TRUE="RangedData") {
  already_checked <- attr(asRangedData, "already_checked")
  if (identical(already_checked, TRUE))
    return(asRangedData)
  if (!isTRUEorFALSE(asRangedData))
    stop("'asRangedData' must be TRUE or FALSE")
  if (asRangedData)
    .Deprecated(msg=.asRangedData_warning_msg(fname, if.FALSE, if.TRUE))
  attr(asRangedData, "already_checked") <- TRUE
  asRangedData
}

