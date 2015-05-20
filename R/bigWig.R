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
  seqlengths <- .Call(BWGFile_seqlengths, path.expand(path(x)))
  Seqinfo(names(seqlengths), seqlengths) # no circularity information
})

setClass("BigWigFileList", contains = "SimpleList",
    prototype = prototype(elementType = "BigWigFile"))

BigWigFileList <- function(path)
{
    new("BigWigFileList", listData = (lapply(path, BigWigFile)))
}

setMethod(path, "BigWigFileList",
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
            object <- as(object, "GRanges")
            callGeneric()
          })

setMethod("export", c("GenomicRanges", "BigWigFile"),
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
                                       seqlengths, compress, con)))
          })

setMethod("export", c("List", "BigWigFile"),
          function(object, con, format, compress = TRUE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            con <- path.expand(path(con))
            if (!isTRUEorFALSE(compress))
              stop("'compress' must be TRUE or FALSE")
            seqlengths <- elementLengths(object)
            sectionPtr <- NULL # keep adding to the same linked list
            on.exit(.Call(BWGSectionList_cleanup, sectionPtr))
            writer <- BigWigWriter(object)
            for(chr in names(object)) {
              sectionPtr <- writer(chr, sectionPtr)
            }
            invisible(BigWigFile(.Call(BWGSectionList_write, sectionPtr,
                                       seqlengths, compress, con)))
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


### FIXME: For chr20, WGS coverage, coercion RangedData->GRanges takes
### 2X as long as reading. The most common use case is reading long
### vectors as Rles. Constructing an RleList is complicated by
### 'which'. As long as 'which' spans chromosomes, we can quickly
### coerce the GRanges (or faster yet, a RangedData) to an
### RleList. Consumes memory, but time cost is minimal, when 'which'
### is trivial (an entire genome or chromosome).

### Consider the use cases:
### - Extract whole genome or whole chromosome coverage as Rle (ChIP-seq)
###   - Use rdToRle() fast path
### - Extract single position coverage for millions of variants
###   - Slow to query by position, so extract whole genome/chromosome (above),
###     then use findRun trick.
### - Extract single position coverage for a hundred hot-spot variants
###   - rdToRle() should be fast enough, followed by findRun trick.
### - Extract coverage for one/all genes and summarize (RNA-seq?)
###   - probably want summarize,BigWigFile


setMethod("import", "BigWigFile",
          function(con, format, text, selection = BigWigSelection(which, ...),
                   which = con, asRangedData = FALSE, asRle = FALSE,
                   as = c("GRanges", "RleList", "NumericList"), ...)
          {
            if (asRangedData) {
              stop("'asRangedData' argument is defunct")
            }
            if (asRle) {
              warning("'asRle' argument has been deprecated, ",
                      "use 'as=\"RleList\"' instead")
              as <- "RleList"
            }
            if (!missing(format))
              checkArgFormat(con, format)
            as <- match.arg(as)

            if (as != "NumericList")
                asRangedData <- normarg_asRangedData(asRangedData, "import")
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
            if (as != "NumericList") {
              which <- as(which, "NormalIRangesList")
            }
            rd <- .Call(BWGFile_query, path.expand(path(con)),
                        as.list(which),
                        identical(colnames(selection), "score"), 
                        as == "NumericList")
            if (as == "NumericList") {
                rd <- as(rd, "NumericList")
                names(rd) <- rep(names(which), elementLengths(which))
                metadata(rd) <- list(ranges = as(which, "GRanges"))
                rd
            } else {
              seqinfo(rd) <- si
              if (as == "RleList") {
                rdToRle(rd)
              } else if (as == "GRanges") {
                strand(rd) <- "*"
                rd <- as(rd, "GRanges")
              } else rd
            }
           })

rdToRle <- function(x) {
  RleList(mapply(function(r, v, sl) {
    coverage(r, width=sl, weight=v$score)
  }, ranges(x), values(x), seqlengths(x), SIMPLIFY=FALSE), compress=FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summary
###

setMethod("summary", "BigWigFile",
          function(object, which = as(seqinfo(object), "GenomicRanges"),
                   size = 1L, type = c("mean", "min", "max", "coverage", "sd"),
                   defaultValue = NA_real_, asRle = FALSE, ...)
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
    dest <- path.expand(dest)
    ans <- .Call(BWGFile_fromWIG, x, seqlengths, dest)
    invisible(BigWigFile(ans))
  }
