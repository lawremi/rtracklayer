### =========================================================================
### BigWig support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BigWigFile", contains = "BiocFile")
setClass("BWFile", contains = "BigWigFile")

BigWigFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("BigWigFile", resource = path)
}
BWFile <- BigWigFile

setMethod("seqinfo", "BigWigFile", function(x) {
  seqlengths <- .Call(BWGFile_seqlengths, expandPath(path(x)))
  Seqinfo(names(seqlengths), seqlengths) # no circularity information
})

setClass("BigWigFileList", contains = "BiocFileList",
    prototype = prototype(elementType = "BigWigFile"))

BigWigFileList <- function(path)
{
    new("BigWigFileList", listData = (lapply(path, BigWigFile)))
}

setMethod("path", "BigWigFileList",
    function(object, ...)
{
    sapply(as.list(object), path)
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

BigWigSelection <- function(ranges=IRangesList(), colnames = "score") {
  if (!is.character(colnames) ||
      (length(colnames) && !identical(colnames, "score")))
    stop("'score' is the only valid column for BigWig")
  if (is.character(ranges))
    new("BigWigSelection", GenomicSelection(ranges, colnames = colnames))
  else {
    if (is(ranges, "BigWigFile"))
      ranges <- seqinfo(ranges)
    new("BigWigSelection", ranges = as(ranges, "IntegerRangesList"),
        colnames = colnames)
  }
}

setAs("IntegerRangesList", "BigWigSelection", function(from) {
  new("BigWigSelection", as(from, "RangedSelection"), colnames = "score")
})

setAs("GenomicRanges", "BigWigSelection", function(from) {
  as(as(from, "IntegerRangesList"), "BigWigSelection")
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
            object <- as(object, "GRanges")
            callGeneric()
          })

setMethod("export", c("GenomicRanges", "BigWigFile"),
          function(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                                  "bedGraph"),
                   compress = TRUE, fixedSummaries = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            con <- path.expand(path(con))
            object <- sortBySeqnameAndStart(object)
            score <- score(object)
            if (isValidScore(score))
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
              if (length(chromData) == 0L)
                return()
              if (any(tail(start(chromData), -1) <= head(end(chromData), -1)))
                stop("BigWig ranges cannot overlap")
              sectionPtr <<- .Call(BWGSectionList_add, sectionPtr,
                                   as.vector(seqnames(chromData)[1]),
                                   as(ranges(chromData), "IRanges"),
                                   as.numeric(score(chromData)), dataFormat)
            }
            dataFormat <- match.arg(dataFormat)
            if (dataFormat == "auto")
              format <- chooseGraphType(object)
            else format <- dataFormat
            on.exit(.Call(BWGSectionList_cleanup, sectionPtr))
            if (format == "bedGraph")
              lapply(split(object, seqnames(object)), .bigWigWriter, con,
                     dataFormat)
            else export.wig(object, con, dataFormat = dataFormat,
                            writer = .bigWigWriter, trackLine = FALSE)
            storage.mode(seqlengths) <- "integer"
            invisible(BigWigFile(.Call(BWGSectionList_write, sectionPtr,
                                       seqlengths, compress, fixedSummaries,
                                       con)))
          })

setMethod("export", c("List", "BigWigFile"),
          function(object, con, format, compress = TRUE, fixedSummaries = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            con <- path.expand(path(con))
            if (!isTRUEorFALSE(compress))
              stop("'compress' must be TRUE or FALSE")
            if (is.null(names(object)))
                stop("'object' must have names")
            seqlengths <- elementNROWS(object)
            sectionPtr <- NULL # keep adding to the same linked list
            on.exit(.Call(BWGSectionList_cleanup, sectionPtr))
            writer <- BigWigWriter(object)
            for(chr in names(object)) {
              sectionPtr <- writer(chr, sectionPtr)
            }
            invisible(BigWigFile(.Call(BWGSectionList_write, sectionPtr,
                                       seqlengths, compress, fixedSummaries,
                                       con)))
          })

setGeneric("BigWigWriter", function(x) standardGeneric("BigWigWriter"))

setMethod("BigWigWriter", "RleList", function(x) {
  function(chr, sectionPtr) {
    .Call(BWGSectionList_add, sectionPtr,
          chr, ranges(x[[chr]]), as.numeric(runValue(x[[chr]])),
          "bedGraph")
  }
})

setMethods("BigWigWriter", list("IntegerList", "NumericList"), function(x) {
  function(chr, sectionPtr) {
    .Call(BWGSectionList_add, sectionPtr,
          chr, NULL, as.numeric(x[[chr]]),
          "fixedStep")
  }
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

### Consider the use cases:
### - Extract whole genome or whole chromosome coverage as Rle (ChIP-seq)
###   - Use coverage() fast path
### - Extract single position coverage for millions of variants
###   - Slow to query by position, so extract whole genome/chromosome (above),
###     then use findRun trick.
### - Extract single position coverage for a hundred hot-spot variants
###   - coverage() should be fast enough, followed by findRun trick.
### - Extract coverage for one/all genes and summarize (RNA-seq?)
###   - probably want summarize,BigWigFile


setMethod("import", "BigWigFile",
          function(con, format, text, selection = BigWigSelection(which, ...),
                   which = con,
                   as = c("GRanges", "RleList", "NumericList"), ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            as <- match.arg(as)
            if (is(which, "GenomicRanges") && as == "NumericList") {
                orig_order <- order(seqnames(which))
            }
            selection <- as(selection, "BigWigSelection")
            validObject(selection)
            si <- seqinfo(con)
            which <- ranges(selection)
            badSpaces <- setdiff(names(which)[lengths(which) > 0L],
                                 seqlevels(si))
            if (length(badSpaces) > 0L)
              warning("'which' contains seqnames not known to BigWig file: ",
                      paste(badSpaces, collapse = ", "))
            which <- which[names(which) %in% seqlevels(si)]
            flatWhich <- unlist(which, use.names = FALSE)
            if (is.null(flatWhich))
              flatWhich <- IRanges()
            which_rl <- split(flatWhich, factor(space(which), seqlevels(si)))
            if (as != "NumericList") {
                which_rl <- as(which, "NormalIRangesList")
            }
            which <- GRanges(which_rl)
            names(which) <- names(unlist(which_rl, use.names=FALSE))
            C_ans <- .Call(BWGFile_query, expandPath(path(con)),
                           as.character(seqnames(which)), ranges(which),
                           identical(colnames(selection), "score"), 
                           as == "NumericList")
            if (as == "NumericList") {
              ans <- as(C_ans, "NumericList")
              names(ans) <- names(which)
              metadata(ans) <- list(ranges = as(which, "GRanges"))
              if (exists("orig_order")) {
                  ans[orig_order] <- ans
              }
              ans
            } else {
              nhits <- C_ans[[3L]]
              gr <- GRanges(rep(seqnames(which), nhits), C_ans[[1L]],
                            seqinfo=si)
              gr$score <- C_ans[[2L]]
              if (as == "RleList") {
                coverage(gr, weight="score")
              } else {
                strand(gr) <- "*"
                gr
              }
            }
           })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summary
###

setMethod("summary", "BigWigFile",
          function(object, which = as(seqinfo(object), "GenomicRanges"),
                   size = 1L, type = c("mean", "min", "max", "coverage", "sd"),
                   defaultValue = NA_real_, asRle = FALSE,
                   as = c("GRangesList", "RleList", "matrix"), ...)
          {
            ### FIXME: could do with "GenomicRanges" here, but
            ### coercions generally only exist for GRanges specifically
            which <- as(which, "GRanges")
            if (!is.numeric(size))
              stop("'size' must be numeric")
            size <- recycleIntegerArg(size, "size", length(which))
            as <- match.arg(as)
            if (any(size > width(which)) && as != "matrix")
              stop("some 'which' are smaller than 'size'; ",
                   "consider passing as='matrix'")
            if (as == "matrix" && (length(size) == 0L || any(size != size[1L])))
                stop("for as='matrix', there must be one unique 'size'")
            type <- match.arg(type)
            if (!missing(asRle))
              warning("argument asRle is deprecated; use as='RleList'")
            if (asRle)
              as <- "RleList"
            if (type == "sd") type <- "std"
            if (!isSingleNumberOrNA(defaultValue))
              stop("'defaultValue' must be a single number or NA")
            summaryList <- .Call(BWGFile_summary, path.expand(path(object)),
                                 as.character(seqnames(which)),
                                 ranges(which), size, type,
                                 as.numeric(defaultValue))
            if (as == "matrix") {
                return(do.call(rbind, summaryList))
            }
            tiles <- tile(which, n = size)
            if (asRle) {
              setNames(RleList(mapply(Rle, summaryList, as.list(width(tiles))),
                               compress=FALSE),
                       names(which))
            } else {
              tiles <- unlist(tiles)
              tiles$score <- unlist(summaryList)
              relist(tiles, summaryList)
            }
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Conversion
###

wigToBigWig <-
  function(x, seqinfo,
           dest = paste(file_path_sans_ext(x, TRUE), "bw", sep = "."),
           clip = FALSE)
  {
    if (!isSingleString(x))
      stop("'x' must be a single string, the path to a WIG file")
    if (!isSingleString(dest))
      stop("'dest' must be a single string, the path to the BigWig output")
    if (!is(seqinfo, "Seqinfo"))
      stop("'seqinfo' must be NULL or a Seqinfo object")
    if (!isTRUEorFALSE(clip))
      stop("'clip' must be TRUE or FALSE")
    seqlengths <- seqlengths(seqinfo)
    if (any(is.na(seqlengths)))
      stop("'seqlengths(seqinfo)' must not contain any 'NA' values")
    x <- path.expand(x)
    dest <- path.expand(dest)
    ans <- .Call(BWGFile_fromWIG, x, clip, seqlengths, dest)
    invisible(BigWigFile(ans))
  }

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

## Remote data are cached locally. Need a way to cleanup.

cleanupBigWigCache <- function(maxDays = 0) {
  stopifnot(isSingleNumber(maxDays))
  dir <- "/tmp/udcCache"
  if (file.exists(dir)) {
      invisible(.Call(R_udcCleanup, maxDays))
  }
}

