### =========================================================================
### UCSC 2bit support
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TwoBitFile class
###

setClass("TwoBitFile", contains = "BiocFile")
setClass("2BitFile", contains = "TwoBitFile")

twoBitPath <- function(path) {
  uri <- .parseURI(path)
  if (!uriIsLocal(uri))
    stop("TwoBit driver handles only local file paths")
  path.expand(uri$path)
}

TwoBitFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("TwoBitFile", resource = twoBitPath(path))
}
`2BitFile` <- TwoBitFile

.seqlengths_TwoBitFile <- function(x) {
    .Call(TwoBitFile_seqlengths, path(x))
}

fastaSeqnames <- function(x) {
    sub(" .*", "", x)
}

setMethod("seqinfo", "TwoBitFile", function(x) {
              seqlengths <- .seqlengths_TwoBitFile(x)
              names(seqlengths) <- fastaSeqnames(names(seqlengths))
              Seqinfo(names(seqlengths), seqlengths)
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

setMethod("export", c("DNAStringSet", "TwoBitFile"),
          function(object, con, format) {
            if (!missing(format))
              checkArgFormat(con, format)  
            seqnames <- names(object)
            if (is.null(seqnames))
              seqnames <- as.character(seq(length(object)))
            freq <- alphabetFrequency(object)
            unsupported.chars <- setdiff(DNA_ALPHABET, c(DNA_BASES, "N"))
            if (any(uniqueLetters(object) %in% unsupported.chars)) {
              stop("One or more strings contain unsupported ambiguity ",
                   "characters.\nStrings can contain only A, C, G, T or N.",
                   "\nSee Biostrings::replaceAmbiguities().")
            }
            if (any(width(object) == 0L)) {
              stop("Empty strings are not yet supported")
            }
            invisible(.TwoBits_export(mapply(.DNAString_to_twoBit, object,
                                             seqnames),
                                      twoBitPath(path(con))))
          })

## Hidden export of a list of twoBit pointers.
## NOT exported (but used in the BSgenome package).
.TwoBits_export <- function(object, con) {
  if (!isSingleString(con))
    stop("'con' must be a single, non-NA string")
  if (!all(sapply(object,
                  function(x) typeof(x) == "externalptr" && is(x, "twoBit"))))
    stop("'object' must be a list of 'twoBit' pointers")
  .Call(TwoBits_write, object, con)
}

## NOT exported (but used in the BSgenome package).
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
            sl <- .seqlengths_TwoBitFile(con)
            sn <- extractROWS(names(sl), match(seqnames(which), seqlevels(con)))
            if (any(is.na(sn))) {
                stop("'seqnames' not in 2bit file: ",
                     paste0("'", unique(seqnames(which)[is.na(sn)]), "'",
                            collapse=", "))
            }
            ans <- .Call(TwoBitFile_read, twoBitPath(path(con)),
                         sn, as(ranges(which), "IRanges"), lkup)
            names(ans) <- names(which)
            ans
          })

setMethod("getSeq", "TwoBitFile",
          function(x, which = as(seqinfo(x), "GenomicRanges")) {
              ans <- import(x, which = which)
              rc <- strand(which) == "-"
              ans[rc] <- reverseComplement(ans[rc])
              ans
          })
