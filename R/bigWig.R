### =========================================================================
### BigWig support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BigWigFile", contains = "RTLFile")
setClass("BWFile", contains = "BigWigFile")

BigWigFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("BigWigFile", resource = path)
}
BWFile <- BigWigFile

setMethod("seqinfo", "BigWigFile", function(x) {
  seqlengths <- .Call(BWGFile_seqlengths, path(x))
  Seqinfo(names(seqlengths), seqlengths) # no circularity information
})

.allowedColNames <- list(bigWig = "score")

.validateColNames <- function(object, format) {
  allowedColNames <- .allowedColNames[[format]]
  invalidColNames <- setdiff(colnames(object), allowedColNames)
  if (length(invalidColNames))
    paste("Column names",
          paste("'", invalidColNames, "'", sep = "", collapse=", "),
          "are invalid for the", format, "format.")
  else NULL
}

setClass("BigWigSelection", prototype = prototype(colnames = "score"),
         contains = "RangedSelection")

setValidity("BigWigSelection",
            function(object) {
              .validateColNames(object, "bigWig")
            })

BigWigSelection <- function(ranges = RangesList(), colnames = "score") {
  if (!is.character(colnames) ||
      (length(colnames) && !identical(colnames, "score")))
    stop("'score' is the only valid column for BigWig")
  if (is.character(ranges))
    new("BigWigSelection", GenomicSelection(ranges, colnames = colnames))
  else {
    if (is(ranges, "BigWigFile"))
      ranges <- seqinfo(ranges)
    new("BigWigSelection", ranges = as(ranges, "RangesList"),
        colnames = colnames)
  }
}

setAs("RangesList", "BigWigSelection", function(from) {
  new("BigWigSelection", as(from, "RangedSelection"), colnames = "score")
})

setAs("GenomicRanges", "BigWigSelection", function(from) {
  as(as(from, "RangesList"), "BigWigSelection")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.bw", function(object, con, ...) standardGeneric("export.bw"))

setMethod("export.bw", "ANY",
          function(object, con, ...)
          {
            export(object, con, "BigWig", ...)
          })

setMethod("export", c("ANY", "BigWigFile"),
          function(object, con, format, ...)
          {
            object <- as(object, "RangedData")
            callGeneric()
          })

setMethod("export", c("RangedData", "BigWigFile"),
          function(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                                  "bedGraph"),
                   compress = TRUE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            con <- path.expand(path(con))
            object <- sortBySeqnameAndStart(object)
            score <- score(object)
            if (!is.numeric(score) || any(is.na(score)))
              stop("The score must be numeric, without any NA's")
            if (!isTRUEorFALSE(compress))
              stop("'compress' must be TRUE or FALSE")
            seqlengths <- seqlengths(object)
            if (any(is.na(seqlengths)))
              stop("Unable to determine seqlengths; either specify ",
                   "'seqlengths' or specify a genome on 'object' that ",
                   "is known to BSgenome or UCSC")
            sectionPtr <- NULL # keep adding to the same linked list
            .bigWigWriter <- function(chromData, con, dataFormat, append) {
              if (any(tail(start(chromData), -1) <= head(end(chromData), -1)))
                stop("BigWig ranges cannot overlap")
              sectionPtr <<- .Call(BWGSectionList_add, sectionPtr,
                                   names(chromData)[1],
                                   as(ranges(chromData)[[1]], "IRanges"),
                                   as.numeric(score(chromData)), dataFormat)
            }
            dataFormat <- match.arg(dataFormat)
            if (dataFormat == "auto")
              format <- chooseGraphType(object)
            else format <- dataFormat
            on.exit(.Call(BWGSectionList_cleanup, sectionPtr))
            if (format == "bedGraph")
              lapply(object, .bigWigWriter, con, dataFormat)
            else export.wig(object, con, dataFormat = dataFormat,
                            writer = .bigWigWriter, trackLine = FALSE)
            storage.mode(seqlengths) <- "integer"
            invisible(BigWigFile(.Call(BWGSectionList_write, sectionPtr,
                                       seqlengths, compress, con)))
          })

setMethod("export", c("RleList", "BigWigFile"),
          function(object, con, format, ...)
          {
            rd <- as(object, "RangedData")
            seqlengths(rd) <- elementLengths(object)
            export(rd, con, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.bw", function(con, ...) standardGeneric("import.bw"))

setMethod("import.bw", "ANY",
          function(con, ...)
          {
            import(con, "BigWig", ...)
          })

setMethod("import", "BigWigFile",
          function(con, format, text, selection = BigWigSelection(which, ...),
                   which = con, asRangedData = TRUE, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            selection <- as(selection, "BigWigSelection")
            validObject(selection)
            si <- seqinfo(con)
            which <- ranges(selection)
            badSpaces <- setdiff(names(which), seqnames(si))
            if (length(badSpaces))
              stop("'which' contains sequence names not known to BigWig file: ",
                   paste(badSpaces, collapse = ", "))
            flatWhich <- unlist(which, use.names = FALSE)
            if (is.null(flatWhich))
              flatWhich <- IRanges()
            which <- split(flatWhich, factor(space(which), seqnames(si)))
            normRanges <- as(which, "NormalIRangesList")
            rd <- .Call(BWGFile_query, path.expand(path(con)),
                        as.list(normRanges),
                        identical(colnames(selection), "score"))
            seqinfo(rd) <- si
            if (!asRangedData) {
              strand(rd) <- "*"
              as(rd, "GRanges")
            }
            else rd
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summary
###

setMethod("summary", "BigWigFile",
          function(object, which = as(seqinfo(object), "GenomicRanges"),
                   size = 1L, type = c("mean", "min", "max", "coverage", "sd"),
                   defaultValue = NA_real_)
          {
            ### FIXME: could do with "GenomicRanges" here, but
            ### coercions generally only exist for GRanges specifically
            which <- as(which, "GRanges")
            if (!is.numeric(size))
              stop("'size' must be numeric")
            size <- recycleIntegerArg(size, "size", length(which))
            type <- match.arg(type)
            if (type == "sd") type <- "std"
            if (!isSingleNumberOrNA(defaultValue))
              stop("'defaultValue' must be a single number or NA")
            summaryList <- .Call(BWGFile_summary, path.expand(path(object)),
                                 as.character(seqnames(which)),
                                 ranges(which), size, type,
                                 as.numeric(defaultValue))
            names(summaryList) <- names(which)
            RleList(summaryList)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Conversion
###

wigToBigWig <-
  function(x, seqinfo,
           dest = paste(file_path_sans_ext(x, TRUE), "bw", sep = "."))
  {
    if (!isSingleString(x))
      stop("'x' must be a single string, the path to a WIG file")
    if (!isSingleString(dest))
      stop("'dest' must be a single string, the path to the BigWig output")
    if (!is(seqinfo, "Seqinfo"))
      stop("'seqinfo' must be NULL or a Seqinfo object")
    seqlengths <- seqlengths(seqinfo)
    if (any(is.na(seqlengths)))
      stop("'seqlengths(seqinfo)' must not contain any 'NA' values")
    x <- path.expand(x)
    ans <- .Call(BWGFile_fromWIG, x, seqlengths, dest)
    invisible(BigWigFile(ans))
  }
