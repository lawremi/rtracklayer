### =========================================================================
### UCSC 2bit support
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TwoBitFile class
###

setClass("TwoBitFile", contains = "RTLFile")
setClass("2BitFile", contains = "TwoBitFile")

twoBitPath <- function(path) {
  uri <- .parseURI(path)
  if (!uriIsLocal(uri))
    stop("TwoBit driver handles only local file paths")
  uri$path
}

TwoBitFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("TwoBitFile", resource = twoBitPath(path))
}
`2BitFile` <- TwoBitFile

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
  export(object, con, "2bit", ...)
})

setMethod("export", c("ANY", "TwoBitFile"), function(object, con, format, ...) {
  object <- as(object, "DNAStringSet")
  callGeneric()
})

setMethod("export", c("BSgenome", "TwoBitFile"),
          function(object, con, format, ...) {
            if (!missing(format))
              checkArgFormat(con, format)
            i <- 0L
            twoBits <- bsapply(new("BSParams", X = object, FUN = function(chr) {
              i <<- i + 1
              .DNAString_to_twoBit(chr, seqnames(object)[i])
            }, ...))
            invisible(.TwoBits_export(as.list(twoBits), twoBitPath(path(con))))
          })

setMethod("export", c("DNAStringSet", "TwoBitFile"),
          function(object, con, format) {
            if (!missing(format))
              checkArgFormat(con, format)  
            seqnames <- names(object)
            if (is.null(seqnames))
              seqnames <- as.character(seq(length(object)))
            invisible(.TwoBits_export(mapply(.DNAString_to_twoBit, object,
                                             seqnames),
                                      twoBitPath(path(con))))
          })

## Hidden export of a list of twoBit pointers
.TwoBits_export <- function(object, con) {
  if (!isSingleString(con))
    stop("'con' must be a single, non-NA string")
  if (!all(sapply(object,
                  function(x) typeof(x) == "externalptr" && is(x, "twoBit"))))
    stop("'object' must be a list of 'twoBit' pointers")
  .Call(TwoBits_write, object, con)
}

.DNAString_to_twoBit <- function(object, seqname) {
  if (!isSingleString(seqname))
    stop("'seqname' must be a single, non-NA string")
  if (!is(object, "DNAString") && !is(object, "MaskedDNAString"))
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

setMethod("import.2bit", "ANY",
          function(con, ...)
          {
            import(con, "2bit", ...)
          })

setMethod("import", "TwoBitFile",
          function(con, format, text, which = as(seqinfo(con), "GenomicRanges"),
                   ...)
          {
            lkup <- get_seqtype_conversion_lookup("B", "DNA")
            ans <- .Call(TwoBitFile_read, as.character(path(con)),
                         as.character(seqnames(which)),
                         as(ranges(which), "IRanges"), lkup)
            names(ans) <- names(which)
            ans
          })
