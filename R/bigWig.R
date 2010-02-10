### UCSC bigWig format 

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

setClass("BigWigSelection", contains = "RangedSelection")

setValidity("BigWigSelection",
            function(object) {
              .validateColNames(object, "bigWig")
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
                   seqlengths, compress, genome, ...)
          {
            rd <- as(object, "RangedData")
            if (!is.null(genome))
              genome(rd) <- genome
            export.bw(rd, con, dataFormat, seqlengths, compress, ...)
          })

setMethod("export.bw", c("RangedData", "character"),
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep",
                                  "bedGraph"),
                   seqlengths, compress)
          {
            score <- score(object)
            if (!is.numeric(score) || any(is.na(score)))
              stop("The score must be numeric, without any NA's")
            if (!IRanges:::isTRUEorFALSE(compress))
              stop("'compress' must be TRUE or FALSE")
            if (is.null(seqlengths) && !is.null(genome(object)))
              seqlengths <- seqlengths(.genomeForID(genome(object)))
            if (is.null(seqlengths))
              stop("Unable to determine seqlengths; either specify ",
                   "'seqlengths' or specify a genome on 'object' that ",
                   "is known to BSgenome")
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

setGeneric("import.bw",
           function(con, selection = GenomicSelection(...), ...)
           standardGeneric("import.bw"))


setMethod("import.bw", "character",
          function(con, selection = GenomicSelection(...), ...)
          {
            if (!IRanges:::isSingleString(con))
              stop("'con' must be a single string, specifying a path")
            selection <- try(as(selection, "RangedSelection"))
            if (is.character(selection))
              stop("'selection' must be coercible to RangedSelection")
            normRanges <- as(ranges(selection), "NormalIRangesList")
            rd <- .Call(BWGFile_query, con, as.list(normRanges),
                        colnames(selection))
            ## Unfortunately, the bigWig query API is such that we can
            ## end up with multiple hits.
            if (any(width(rd) > 2)) {
              hits <- queryHits(findOverlaps(ranges(rd), normRanges))
              rd <- rd[!duplicated(hits),]
            }
            rd
          })
