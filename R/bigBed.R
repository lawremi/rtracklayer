### =========================================================================
### BigBed support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BigBedFile", contains = "RTLFile")
setClass("BBFile", contains = "BigBedFile")

BigBedFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("BigBedFile", resource = path)
}
BBFile <- BigBedFile

setMethod("seqinfo", "BigBedFile", function(x) {
  seqlengths <- .Call(BBDFile_seqlengths, expandPath(path(x)))
  Seqinfo(names(seqlengths), seqlengths)
})

.defaultColNames <- c("name", "score", "thick", "itemRgb", "blocks")

setClass("BigBedSelection", prototype = prototype(colnames = .defaultColNames),
         contains = "RangedSelection")

BigBedSelection <- function(ranges=IRangesList(), colnames = .defaultColNames) {
  if (is.character(ranges))
    new("BigBedSelection", GenomicSelection(ranges, colnames = colnames))
  else {
    if (is(ranges, "BigBedFile"))
      ranges <- seqinfo(ranges)
    new("BigBedSelection", ranges = as(ranges, "IntegerRangesList"),
     colnames = colnames)
  }
}
