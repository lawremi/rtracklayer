### =========================================================================
### BED (Browser Extended Display) support (including bedGraph and BED15)
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BEDFile", contains = "RTLFile")
BEDFile <- function(resource) {
  new("BEDFile", resource = resource)
}

setClass("BEDGraphFile", contains = "BEDFile")
BEDGraphFile <- function(resource) {
  new("BEDGraphFile", resource = resource)
}

setClass("BED15File", contains = "BEDFile")
BED15File <- function(resource) {
  new("BED15File", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.bed",
           function(object, con, ...) standardGeneric("export.bed"))

setMethod("export.bed", "ANY",
          function(object, con, ...) {
            export(object, con, "bed", ...)
          })

setMethod("export", c("ANY", "BEDFile"),
          function(object, con, format, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            cl <- class(object)
            if (hasMethod("asBED", class(object)))
              object <- asBED(object)
            track <- try(as(object, "RangedData"), silent = TRUE)
            if (class(track) == "try-error") {
              track <- try(as(object, "RangedDataList"), silent = TRUE)
              if (class(track) == "try-error")
                stop("cannot export object of class '", cl, "'")
            }
            export(track, con, ...)
          })

setMethod("export", c("RangedData", "BEDFile"),
          function(object, con, format, append = FALSE, index = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            name <- strand <- thickStart <- thickEnd <- color <- NULL
            blockCount <- blockSizes <- blockStarts <- NULL
            if (index)
              object <- sortBySeqnameAndStart(object)
            df <- data.frame(seqnames(object), start(object) - 1L, end(object))
            score <- score(object)
            if (!is.null(score)) {
              if (!is.numeric(score) || any(is.na(score)))
                stop("Scores must be non-NA numeric values")
            }
            if (is(con, "BEDGraphFile")) {
              if (is.null(score)) ## bedGraph requires score
                score <- 0
              df$score <- score
            } else {
              toCSV <- function(x) {
                if (is(x, "IntegerList")) {
                  x <- unlist(lapply(x, paste, collapse = ","), use.names=FALSE)
                } else if (!is.character(x) && !is.null(x))
                  stop("Could not convert block coordinates to CSV")
                x
              }
              blockSizes <- blockStarts <- NULL
              if (!is.null(object$blocks)) {
                blockSizes <- toCSV(width(object$blocks))
                blockStarts <- toCSV(start(object$blocks) - 1L)
              }
              if (is(con, "BED15File") && is.null(blockSizes))
                blockStarts <- blockSizes <- "" # bed15 must have all 15 cols
              if (!is.null(blockSizes) || !is.null(blockStarts)) {
                if (is.null(blockSizes))
                  stop("'blockStarts' specified without 'blockSizes'")
                if (is.null(blockStarts))
                  stop("'blockSizes' specified without 'blockStarts'")
                lastBlock <- function(x)
                  as.integer(sub(".*,", "", sub(",$", "", x)))
                lastSize <- lastBlock(blockSizes)
                lastStart <- lastBlock(blockStarts)
                if (any(df[[2]] + lastSize + lastStart != df[[3]]) ||
                    any(sub(",.*", "", blockStarts) != "0"))
                  stop("blocks must span entire feature")
                blockCount <- elementLengths(object$blocks)
                if (!is.null(object$blockCount))
                  if (!identical(blockCount, as.integer(object$blockCount)))
                    stop("incorrect block counts given block sizes")
              }
              if (is.null(color))
                color <- object$itemRgb
              if (is.null(color) && !is.null(blockCount))
                color <- "0" ## blocks require color
              else if (!is.null(color)) {
                nacol <- is.na(color)
                colmat <- col2rgb(color)
                color <- paste(colmat[1,], colmat[2,], colmat[3,], sep = ",")
                color[nacol] <- "0"
              }
              thickStart <- object$thickStart
              thickEnd <- object$thickEnd
              if (!is.null(object$thick)) {
                thickStart <- start(object$thick)
                thickEnd <- end(object$thick)
              }
              ## color requires thick ranges
              if (is.null(thickStart) && !is.null(color)) {
                thickStart <- start(object)
                thickEnd <- end(object)
              }
              strand <- object$strand
              if (!is.null(thickStart) && is.null(strand)) {
                strand <- rep(NA, nrow(object))
              }
              if (!is.null(strand) && is.null(score))
                score <- 0
              name <- object$name
              if (is.null(name))
                name <- rownames(object)
              if (!is.null(score) && is.null(name))
                name <- rep(NA, nrow(object))
              df$name <- name
              df$score <- score
              df$strand <- strand
              df$thickStart <-
                if (!is.null(thickStart)) thickStart - 1L else NULL
              df$thickEnd <- thickEnd
              df$itemRgb <- color
              df$blockCount <- blockCount
              df$blockSizes <- blockSizes
              df$blockStarts <- blockStarts
              if (is(con, "BED15File")) {
                df$expCount <- object$expCount
                df$expIds <- object$expIds
                df$expScores <- object$expScores
              }
            }
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            file <- con
            con <- connection(con, if (append) "a" else "w")
            on.exit(release(con))
            write.table(df, con, sep = "\t", col.names = FALSE,
                        row.names = FALSE, quote = FALSE, na = ".")
            release(con)
            if (index)
              invisible(indexTrack(file))
            else invisible(file)
          })

setMethod("export", c("UCSCData", "BEDFile"),
          function(object, con, format, trackLine = TRUE, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (trackLine) {
              export.ucsc(object, con, ...)
            } else {
              callNextMethod()
            }
            invisible(con)
          })

setMethod("export", c("RangedDataList", "BEDFile"),
          .export_RangedDataList_RTLFile)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.bed", function(con, ...) standardGeneric("import.bed"))

setMethod("import.bed", "ANY",
          function(con, ...)
          {
            import(con, format = "bed", ...)
          })

scanTrackLine <- function(con) {
  con <- connectionForResource(con, "r")
  line <- "#"
  while(length(grep("^ *#", line))) # skip initial comments
    line <- readLines(con, 1, warn = FALSE)
  if (length(grep("^track", line)) > 0)
    line
  else {
    pushBack(line, con)
    NULL
  }
}

setMethod("import", "BEDFile",
          function(con, format, text, trackLine = TRUE,
                   genome = NA, asRangedData = TRUE, colnames = NULL,
                   which = NULL, seqinfo = NULL, extraCols = character())
          {
            if (!missing(format))
              checkArgFormat(con, format)
            file <- con
            con <- queryForConnection(con, which)
            if (attr(con, "usedWhich"))
              which <- NULL
            if (is.null(seqinfo))
              seqinfo <- attr(con, "seqinfo")
            ## check for a track line
            line <- scanTrackLine(con)
            if (!is.null(line) && trackLine) {
              pushBack(line, con)
              return(import.ucsc(initialize(file, resource = con), drop = TRUE,
                                 trackLine = FALSE, genome = genome,
                                 asRangedData = asRangedData,
                                 colnames = colnames,
                                 which = which, seqinfo = seqinfo))
            }
            if (is(file, "BEDGraphFile")) {
              bedClasses <- c("character", "integer", "integer", "numeric")
              bedNames <- c("chrom", "start", "end", "score")
            } else {
              bedNames <- c("chrom", "start", "end", "name",
                            "score", "strand", "thickStart",
                            "thickEnd", "itemRgb", "blockCount",
                            "blockSizes", "blockStarts")
              bedClasses <- c("character", "integer", "integer", "character",
                              "numeric", "character", "integer", "integer",
                              "character", "integer", "character", "character")
            }
            if (is(file, "BED15File"))
              bedNames <- c(bedNames, "expCount", "expIds", "expScores")
            normArgColnames <- function(validNames) {
              if (is.null(colnames))
                colnames <- validNames
              else {
                colnames <- unique(c(head(bedNames, 3), as.character(colnames)))
                if ("thick" %in% colnames)
                  colnames <- c(setdiff(colnames, "thick"), "thickStart",
                                "thickEnd")
                if ("blocks" %in% colnames)
                  colnames <- c(setdiff(colnames, "blocks"), "blockStarts",
                                "blockSizes", "blockCount")
                missingCols <- setdiff(colnames, validNames)
                if (length(missingCols))
                  stop("Requested column(s) ",
                       paste("'", missingCols, "'", sep = "", collapse = ", "),
                       " are not valid columns or were not found in the file")
              }
              colnames
            }
            ## read a single line to get ncols up-front,
            ## and thus specify all col classes
            ## FIXME: reading in 'as.is' to save memory,
            line <- readLines(con, 1, warn=FALSE)
            ## UCSC seems to use '#' at beginning to indicate comment.
            while(length(line) &&
                  (!nzchar(line) || substring(line, 1, 1) == "#"))
            {
              line <- readLines(con, 1, warn=FALSE)
            }
            if (length(line)) {
              `tail<-` <- function(x, n, value)
                if (n != 0) c(head(x, -n), value) else x
              pushBack(line, con)
              colsInFile <- seq_len(length(strsplit(line, "[\t ]")[[1]]))
              presentNames <- bedNames[colsInFile]
              tail(presentNames, length(extraCols)) <- names(extraCols)
              presentClasses <- bedClasses[colsInFile]
              tail(presentClasses, length(extraCols)) <- unname(extraCols)
              colnames <- normArgColnames(presentNames)
              bedClasses <- ifelse(presentNames %in% colnames,
                                   presentClasses, "NULL")
              bed <- DataFrame(read.table(con, colClasses = bedClasses,
                                          as.is = TRUE, na.strings = ".",
                                          comment.char = ""))
            } else {
              if (is.null(colnames))
                colnames <- character()
              else colnames <- normArgColnames(bedNames)
              keepCols <- bedNames %in% colnames
              bed <- DataFrame(as.list(sapply(bedClasses[keepCols], vector)))
            }
            colnames(bed) <- bedNames[bedNames %in% colnames]
            bed <- bed[substring(bed$chrom, 1, 1) != "#",]
            if (!is.null(bed$thickStart)) {
              thickEnd <- bed$thickEnd
              if (is.null(thickEnd))
                thickEnd <- bed$end
              bed$thick <- IRanges(bed$thickStart + 1L, thickEnd)
              bed$thickStart <- bed$thickEnd <- NULL
            }
            color <- bed$itemRgb
            if (is.character(color)) { # could be NULL
              spec <- color != "0"
              cols <- unlist(strsplit(color[spec], ",", fixed=TRUE),
                             use.names=FALSE)
              cols <- matrix(as.integer(cols), 3)
              color <- rep(NA, nrow(bed))
              color[spec] <- rgb(cols[1,], cols[2,], cols[3,], max = 255)
              bed$itemRgb <- color              
            }
            fromCSV <- function(b) {
              as.integer(unlist(strsplit(b, ",", fixed = TRUE)))
            }
            if (!is.null(bed$blockStarts)) {
              blocks <- seqsplit(IRanges(fromCSV(bed$blockStarts) + 1L,
                                         width = fromCSV(bed$blockSizes)),
                                 togroup(PartitioningByWidth(bed$blockCount)))
              names(blocks) <- NULL
              bed$blockStarts <- bed$blockSizes <- bed$blockCount <- NULL
              bed$blocks <- blocks
            }
            GenomicData(IRanges(bed$start + 1L, bed$end),
                        bed[-(1:3)],
                        chrom = bed$chrom, genome = genome,
                        seqinfo = seqinfo,
                        asRangedData = asRangedData, which = which)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BED15 (Microarray) Support
###

setMethod("import", "BED15File",
          function(con, format, text, trackLine = NULL, genome = NA,
                   asRangedData = TRUE, which = NULL)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!isTRUEorFALSE(asRangedData))
              stop("'asRangedData' must be TRUE or FALSE")
            if (is.null(trackLine))
              return(import.ucsc(con, TRUE, genome = genome,
                                 asRangedData = asRangedData, which = which))
            bed <- callNextMethod()
            if (asRangedData) {
              if (!nrow(bed))
                return(bed)
              ids <- strsplit(bed$expIds[1], ",", fixed=TRUE)[[1]]
              expNames <- trackLine@expNames[as.integer(ids) + 1L]
              scores <- unlist(strsplit(bed$expScores, ",", fixed=TRUE),
                               use.names=FALSE)
              scores <- as.numeric(scores)
              scores[scores == -10000] <- NA # stupid UCSC convention
              scores <- split(scores, gl(length(expNames), 1, length(scores)))
              names(scores) <- expNames
              nonExpCols <- setdiff(colnames(bed),
                                    c("expCount", "expScores", "expIds"))
              bed <- bed[,nonExpCols]
              for (samp in names(scores))
                bed[[samp]] <- scores[[samp]]
              bed
            } else {
              if (!length(bed))
                return(bed)
              ids <- strsplit(values(bed)$expIds[1], ",", fixed=TRUE)[[1]]
              expNames <- trackLine@expNames[as.integer(ids) + 1L]
              scores <- unlist(strsplit(values(bed)$expScores, ",", fixed=TRUE),
                               use.names=FALSE)
              scores <- as.numeric(scores)
              scores[scores == -10000] <- NA # stupid UCSC convention
              scores <- split(scores, gl(length(expNames), 1, length(scores)))
              names(scores) <- expNames
              nonExpCols <- setdiff(colnames(values(bed)),
                                    c("expCount", "expScores", "expIds"))
              values(bed) <- values(bed)[,nonExpCols]
              values(bed) <- cbind(values(bed), do.call(DataFrame, scores))
              bed              
            }
          })

setGeneric("import.bed15",
           function(con, ...)
           standardGeneric("import.bed15"),
           signature = "con")

setMethod("import.bed15", "ANY",
          function(con, ...)
          {
            import(con, "bed15", ...)
          })

setGeneric("export.bed15",
           function(object, con, ...) standardGeneric("export.bed15"))

setMethod("export.bed15", "ANY",
          function(object, con, ...) {
            export(object, con, "bed15",  ...)
          })

### FIXME: dispatch will break when 'object' is a UCSCData
### Possible solution: just merge this code with the main BEDFile method?
setMethod("export", c("RangedData", "BED15File"),
          function(object, con, format, expNames = NULL, trackLine = NULL, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (is.null(trackLine)) {
              ## ensure we do not override existing track line parameter
              if (is.null(expNames) && is(object, "UCSCData") &&
                  is(object@trackLine, "Bed15TrackLine"))
                expNames <- object@trackLine@expNames
              return(export.ucsc(object, con, expNames = expNames, ...))
            }
            expNames <- trackLine@expNames
            object$expCount <- rep(length(expNames), nrow(object))
            object$expIds <- rep(paste(seq_along(expNames)-1, collapse=","),
                                 nrow(object))
            scores <- as.list(unlist(values(object[,expNames])))
            scores <- do.call(paste, c(scores, sep = ","))
            scores <- gsub("NA", "-10000", scores, fixed=TRUE)
            object$expScores <- scores
            callNextMethod(object, con, ...)
          })

setClass("Bed15TrackLine",
         representation(expStep = "numeric", expScale = "numeric",
                        expNames = "characterORNULL"),
         prototype(expStep = 0.5, expScale = 3.0), 
         contains = "BasicTrackLine") # not sure which fields work

setAs("Bed15TrackLine", "character",
      function(from)
      {
        str <- as(as(from, "TrackLine"), "character")
        paste(str, " type=array ",
              " expScale=", from@expScale,
              " expStep=", from@expStep,
              " expNames=\"", paste(from@expNames, collapse=","), "\"",
              sep = "")
      })

setAs("character", "Bed15TrackLine",
      function(from)
      {
        line <- new("Bed15TrackLine", as(from, "TrackLine"))
        vals <- ucscParsePairs(from)
        line@expScale <- as.numeric(vals[["expScale"]])
        line@expStep <- as.numeric(vals[["expStep"]])
        line@expNames <- strsplit(vals[["expNames"]], ",", fixed=TRUE)[[1]]
        line
      })

setMethod("fileFormat", "Bed15TrackLine", function(x) "bed15")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bedGraph (formerly subset of WIG) support
###

setGeneric("import.bedGraph",
           function(con, ...) standardGeneric("import.bedGraph"),
           signature = "con")

setMethod("import.bedGraph", "ANY",
          function(con, ...)
          {
            import(con, "bedGraph", ...)
          })

setGeneric("export.bedGraph",
           function(object, con, ...)
           standardGeneric("export.bedGraph"))

setMethod("export.bedGraph", "ANY",
          function(object, con, ...)
          {
            export(object, con, "bedGraph", ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setGeneric("asBED", function(x, ...) standardGeneric("asBED"))

setMethod("asBED", "GRangesList", function(x) {
  x_range <- range(x)
  if (any(elementLengths(x_range) != 1L))
    stop("Empty or multi-strand/seqname elements not supported by BED")
  gr <- unlist(x_range, use.names=FALSE)
  values(gr) <- values(x)
  values(gr)$name <- names(x)
  x_ranges <- ranges(unlist(x, use.names=FALSE))
  ord_start <- order(start(x_ranges))
  x_ranges <- shift(x_ranges, 1L - rep(start(gr), elementLengths(x)))[ord_start]
  values(gr)$blocks <- split(x_ranges, togroup(x)[ord_start])
  gr
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setGeneric("blocks", function(x, ...) standardGeneric("blocks"))

setMethod("blocks", "RangedData",
          function(x)
          {
            blocks(as(x, "GRanges"))
          })

setMethod("blocks", "GenomicRanges",
          function(x)
          {
            block_counts <- elementLengths(values(x)$blocks)
            gr <- GRanges(rep(seqnames(x), block_counts),
                          shift(unlist(values(x)$blocks, use.names = FALSE),
                                rep(start(x), block_counts) - 1L),
                          rep(strand(x), block_counts))
            seqinfo(gr) <- seqinfo(x)
            split(gr, togroup(values(x)$blocks))
          })
