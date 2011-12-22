### UCSC bigWig format 

## NOTE: could use an externalptr here, but we should profile to see
## if that is worth it.
setClass("BigWigFile", contains = "RTLFile")

BigWigFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("BigWigFile", path = path)
}

setMethod("seqinfo", "BigWigFile", function(x) {
  seqlengths <- .Call(BWGFile_seqlengths, x@path)
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

BigWigSelection <- function(ranges = GRanges(), colnames = "score") {
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
  new("BigWigSelection", as(from, "RangedSelection"))
})

setAs("GenomicRanges", "BigWigSelection", function(from) {
  as(as(from, "RangesList"), "BigWigSelection")
})

setGeneric("export.bw",
           function(object, con,
                    dataFormat = c("auto", "variableStep", "fixedStep",
                      "bedGraph"),
                    compress = TRUE, ...)
           standardGeneric("export.bw"))

setMethod("export.bw", "ANY",
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                     "bedGraph"), compress, ...)
          {
            rd <- as(object, "RangedData")
            export.bw(rd, con, dataFormat, compress, ...)
          })

setMethod("export.bw", c("RangedData", "character"),
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                                  "bedGraph"),
                   compress)
          {
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
            .bigWigWriter <- function(chromData, con, dataFormat) {
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
            else export.wigLines(object, con, dataFormat, .bigWigWriter)
            storage.mode(seqlengths) <- "integer"
            invisible(BigWigFile(.Call(BWGSectionList_write, sectionPtr,
                                       seqlengths, compress, con)))
          })

setGeneric("import.bw", function(con, ...) standardGeneric("import.bw"))

setMethod("import.bw", "connection",
          function(con, ...)
          {
            import.bw(summary(con)$description, ...)
          })

setMethod("import.bw", "character",
          function(con, ...)
          {
            import.bw(BigWigFile(con), ...)
          })

setMethod("import.bw", "BigWigFile",
          function(con, selection = BigWigSelection(ranges, ...),
                   ranges = con, asRangedData = TRUE, ...)
          {
            selection <- as(selection, "BigWigSelection")
            validObject(selection)
            normRanges <- as(ranges(selection), "NormalIRangesList")
            rd <- .Call(BWGFile_query, path(con), as.list(normRanges),
                        identical(colnames(selection), "score"))
            ## Unfortunately, the bigWig query API is such that we can
            ## end up with multiple hits.
            if (any(width(rd) > 2)) {
              hits <- queryHits(findOverlaps(ranges(rd), normRanges))
              rd <- rd[!duplicated(hits),]
            }
            if (!asRangedData)
              as(rd, "GRanges")
            else rd
          })

setMethod("summary", "BigWigFile",
          function(object, ranges = as(seqinfo(object), "GenomicRanges"),
                   size = 1L, type = c("mean", "min", "max", "coverage", "sd"),
                   defaultValue = NA_real_)
          {
            ### FIXME: could do with "GenomicRanges" here, but
            ### coercions generally only exist for GRanges specifically
            ranges <- as(ranges, "GRanges")
            if (!is.numeric(size))
              stop("'size' must be numeric")
            size <- recycleIntegerArg(size, "size", length(ranges))
            type <- match.arg(type)
            if (type == "sd") type <- "std"
            if (!isSingleNumberOrNA(defaultValue))
              stop("'defaultValue' must be a single number or NA")
            summaryList <- .Call(BWGFile_summary, path(object),
                                 as.character(seqnames(ranges)),
                                 ranges(ranges), size, type,
                                 as.numeric(defaultValue))
            names(summaryList) <- names(ranges)
            RleList(summaryList)
          })

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
    ans <- .Call(BWGFile_fromWIG, x, seqlengths, dest)
    invisible(BigWigFile(ans))
  }
