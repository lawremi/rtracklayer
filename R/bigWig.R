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
                    seqlengths = NULL, compress = TRUE, ...)
           standardGeneric("export.bw"))

setMethod("export.bw", "ANY",
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                     "bedGraph"),
                   seqlengths, compress, genome = NULL, ...)
          {
            rd <- as(object, "RangedData")
            export.bw(rd, con, dataFormat, seqlengths, compress,
                      genome = genome, ...)
          })

setMethod("export.bw", c("RangedData", "character"),
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                                  "bedGraph"),
                   seqlengths, compress, genome = NULL)
          {
            score <- score(object)
            if (!is.numeric(score) || any(is.na(score)))
              stop("The score must be numeric, without any NA's")
            if (!isTRUEorFALSE(compress))
              stop("'compress' must be TRUE or FALSE")
            if (!is.null(genome))
              genome(object) <- genome
            if (is.null(seqlengths) && !is.null(genome(object))) {
              si <- seqinfoForGenome(singleGenome(genome(object)))
              if (is.null(si))
                stop("Unable to determine seqlengths; either specify ",
                     "'seqlengths' or specify a genome on 'object' that ",
                     "is known to BSgenome or UCSC")
              seqlengths <- seqlengths(si)
            }
            if (!is.numeric(seqlengths) ||
                !all(names(object) %in% names(seqlengths)))
              stop("seqlengths must be numeric and indicate a length for ",
                   "each sequence in 'object'")
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
            .Call(BWGSectionList_write, sectionPtr, seqlengths,
                  compress, con)
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
                   ranges = con, ...)
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
            rd
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
