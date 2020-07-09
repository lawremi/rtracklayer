### =========================================================================
### BigBed support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BigBedFile", contains = "RTLFile")

BigBedFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("BigBedFile", resource = path)
}

setMethod("seqinfo", "BigBedFile", function(x) {
  seqlengths <- .Call(BBDFile_seqlengths, expandPath(path(x)))
  Seqinfo(names(seqlengths), seqlengths)
})
