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
