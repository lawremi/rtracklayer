### =========================================================================
### UCSC 2bit support
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TwoBitFile class
###

setClass("TwoBitFile", contains = "RTLFile")

TwoBitFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("TwoBitFile", path = path)
}

setMethod("seqinfo", "TwoBitFile", function(x) {
  seqlengths <- .Call(TwoBitFile_seqlengths, path(x))
  Seqinfo(names(seqlengths), seqlengths) # no circularity or genome information
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.2bit",
           function(object, con, ...) standardGeneric("export.2bit"))

setMethod("export.2bit", "ANY", function(object, con, ...) {
  export.2bit(as(object, "DNAStringSet"), con, ...)
})

setMethod("export.2bit", c("BSgenome", "character"),
          function(object, con, exclude = character(0), maskList = logical(0)) {
            i <- 0L
            twoBits <- bsapply(new("BSParams", X = object, FUN = function(chr) {
              i <<- i + 1
              .DNAString_to_twoBit(chr, seqnames(object)[i])
            }, exclude = exclude, maskList = maskList))
            .TwoBits_export(as.list(twoBits), con)
          })

setMethod("export.2bit", c("DNAStringSet", "character"), function(object, con) {
  seqnames <- names(object)
  if (is.null(seqnames))
    seqnames <- as.character(seq(length(object)))
  .TwoBits_export(mapply(.DNAString_to_twoBit, object, seqnames), con)
})

## Hidden export of a list of twoBit pointers
.TwoBits_export <- function(object, con) {
  if (!IRanges:::isSingleString(con))
    stop("'con' must be a single, non-NA string")
  if (!all(sapply(object,
                  function(x) typeof(x) == "externalptr" && is(x, "twoBit"))))
    stop("'object' must be a list of 'twoBit' pointers")
  .Call(TwoBits_write, object, con)
}

.DNAString_to_twoBit <- function(object, seqname) {
  if (!IRanges:::isSingleString(seqname))
    stop("'seqname' must be a single, non-NA string")
  if (!is(object, "DNAString"))
    stop("'object' must be a DNAString")
  object_masks <- masks(object)
  if (!is.null(object_masks)) {
    object_mask <- collapse(object_masks)[[1]]
    active(masks(object)) <- FALSE # so we get the real sequence
  } else object_mask <- IRanges()
  .Call(DNAString_to_twoBit, as.character(object), object_mask, seqname)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.2bit", function(con, ...) standardGeneric("import.2bit"))

setMethod("import.2bit", "connection",
          function(con, ...)
          {
            import.2bit(summary(con)$description, ...)
          })

setMethod("import.2bit", "character",
          function(con, ...)
          {
            import.2bit(TwoBitFile(con), ...)
          })

setMethod("import.2bit", "TwoBitFile",
          function(con, which = as(seqinfo(con), "GenomicRanges"), ...)
          {
            lkup <- get_xsbasetypes_conversion_lookup("B", "DNA")
            ans <- .Call(TwoBitFile_read, as.character(path(con)),
                         as.character(seqnames(which)),
                         as(ranges(which), "IRanges"), lkup)
            names(ans) <- names(which)
            ans
          })
