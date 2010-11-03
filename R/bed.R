# Import/export of Browser Extended Display (BED) data

setGeneric("export.bed",
           function(object, con, variant = c("base", "bedGraph", "bed15"),
                    color = NULL, append = FALSE, ...)
           standardGeneric("export.bed"))


setMethod("export.bed", "ANY",
          function(object, con, variant = c("base", "bedGraph", "bed15"), color,
                   append)
          {
            cl <- class(object)
            track <- try(as(object, "RangedData"), silent = TRUE)
            if (class(track) == "try-error") {
              track <- try(as(object, "RangedDataList"), silent = TRUE)
              if (class(track) == "try-error")
                stop("cannot export object of class '", cl, "'")
            }
            export.bed(track, con=con, variant=variant, color=color,
                       append=append)
          })

setMethod("export.bed", c("RangedData", "characterORconnection"),
          function(object, con, variant = c("base", "bedGraph", "bed15"), color,
                   append)
          {
            variant <- match.arg(variant)
            name <- strand <- thickStart <- thickEnd <- color <- NULL
            blockCount <- blockSizes <- blockStarts <- NULL
            df <- data.frame(chrom(object), start(object)-1, end(object))
            score <- score(object)
            if (!is.null(score)) {
              if (!is.numeric(score) || any(is.na(score)))
                stop("Scores must be non-NA numeric values")
            }
            if (variant == "bedGraph") {
              if (is.null(score)) ## bedGraph requires score
                score <- 0
              df$score <- score
            } else {
              blockSizes <- object$blockSizes
              blockStarts <- object$blockStarts
              if (variant == "bed15" && is.null(blockSizes))
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
                blockCount <- length(strsplit(blockSizes, ",")[[1]])
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
              thickEnd <- object$thickEnd ## color requires thick ranges
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
              df$thickStart <- thickStart
              df$thickEnd <- thickEnd
              df$itemRgb <- color
              df$blockCount <- blockCount
              df$blockSizes <- blockSizes
              df$blockStarts <- blockStarts
              if (variant == "bed15") {
                df$expCount <- object$expCount
                df$expIds <- object$expIds
                df$expScores <- object$expScores
              }
            }
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            write.table(df, con, sep = "\t", col.names = FALSE,
                        row.names = FALSE, quote = FALSE, na = ".",
                        append = append)
          })

setMethod("export.bed", "UCSCData",
          function(object, con, variant = c("base", "bedGraph", "bed15"), color,
                   append, trackLine = TRUE, ...)
          {
            variant <- match.arg(variant)
            if (variant == "base" && trackLine) {
              export.ucsc(object, con, "bed", append, variant, color, ...)
            } else {
              callNextMethod(object, con, variant, color, append)
            }
          })

setMethod("export.bed", "RangedDataList",
          function(object, con, variant = c("base", "bedGraph", "bed15"), color,
                   append, ...)
          {
            export.ucsc(object, con, "bed", append, variant, color, ...)
          })

setGeneric("import.bed",
           function(con, variant = c("base", "bedGraph", "bed15"),
                    trackLine = TRUE, genome = "hg18", asRangedData = TRUE,
                    ...)
           standardGeneric("import.bed"))

setMethod("import.bed", "character",
          function(con, variant = c("base", "bedGraph", "bed15"), trackLine,
                   genome, asRangedData = TRUE)
          {
            import(file(con), format = "bed", variant = variant,
                   trackLine = trackLine, genome = genome,
                   asRangedData = asRangedData)
          })

setMethod("import.bed", "connection",
          function(con, variant = c("base", "bedGraph", "bed15"), trackLine,
                   genome, asRangedData = TRUE)
          {
            variant <- match.arg(variant)
            if (variant == "base" && trackLine) {
              ## check for a track line
              line <- "#"
              while(length(grep("^ *#", line))) # skip initial comments
                line <- readLines(con, 1, warn = FALSE)
              pushBack(line, con)
              if (length(grep("^track", line)) > 0)
                return(import.ucsc(con, subformat = "bed", drop = TRUE,
                                   trackLine = FALSE, genome = genome,
                                   asRangedData = asRangedData))
            }
            if (variant == "bedGraph") {
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
            if (variant == "bed15")
              bedNames <- c(bedNames, "expCount", "expIds", "expScores")
            ## read a single line to get ncols up-front,
            ## and thus specify all col classes
            ## FIXME: reading in 'as.is' to save memory,
            if (length(line <- readLines(con, 1, warn=FALSE))) {
              pushBack(line, con)
              bedClasses <-
                bedClasses[seq(length(strsplit(line, "[\t ]")[[1]]))]
              bed <- DataFrame(read.table(con, colClasses = bedClasses,
                                          as.is = TRUE))
            } else bed <- DataFrame(as.list(sapply(bedClasses, vector)))
            colnames(bed) <- bedNames[seq_len(ncol(bed))]
            if (variant != "bedGraph") { ## don't know column #, coerce here
              bed$start <- as.integer(bed$start)
              bed$end <- as.integer(bed$end)
            } ## BED is 0-start, so add 1 to start
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
            GenomicData(IRanges(bed$start + 1L, bed$end),
                        bed[,tail(colnames(bed), -3),drop=FALSE],
                        chrom = bed$chrom, genome = genome,
                        asRangedData = asRangedData)
          })

setGeneric("blocks", function(x, ...) standardGeneric("blocks"))

setMethod("blocks", "RangedData",
          function(x)
          {
            if (is.null(x$blockStarts) || is.null(x$blockSizes))
              stop("'x' must have 'blockStarts' and 'blockSizes' columns")
            starts <- unlist(strsplit(as.character(x$blockStarts), ","))
            sizes <- unlist(strsplit(as.character(x$blockSizes), ","))
            GRanges(rep(chrom(x), x$blockCount),
                    IRanges(as.integer(starts) + rep(start(x), x$blockCount),
                            width = as.integer(sizes)),
                    tx = rep(x$name, x$blockCount))
          })

setGeneric("import.bed15Lines",
           function(con, trackLine, genome = "hg18", asRangedData = TRUE, ...)
           standardGeneric("import.bed15Lines"))

setMethod("import.bed15Lines", "ANY",
          function(con, trackLine, genome, asRangedData = TRUE)
          {
            if (!IRanges:::isTRUEorFALSE(asRangedData))
              stop("'asRangedData' must be TRUE or FALSE")
            bed <- import.bed(con, "bed15", genome = genome,
                              asRangedData = asRangedData)
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
           function(con, genome = "hg18", asRangedData = TRUE, ...)
           standardGeneric("import.bed15"))

setMethod("import.bed15", "ANY",
          function(con, genome, asRangedData = TRUE)
          {
            import.ucsc(con, "bed15", TRUE, genome = genome,
                        asRangedData = asRangedData)
          })

setGeneric("export.bed15Lines",
           function(object, con, trackLine, ...)
           standardGeneric("export.bed15Lines"))

setMethod("export.bed15Lines", "RangedData",
          function(object, con, trackLine, ...)
          {
            expNames <- trackLine@expNames
            object$expCount <- rep(length(expNames), nrow(object))
            object$expIds <- rep(paste(seq_along(expNames)-1, collapse=","),
                                 nrow(object))
            scores <- as.list(unlist(values(object[,expNames])))
            scores <- do.call("paste", c(scores, sep = ","))
            scores <- gsub("NA", "-10000", scores, fixed=TRUE)
            object$expScores <- scores
            export.bed(object, con, "bed15", ...)
          })

setGeneric("export.bed15",
           function(object, con, expNames = NULL, ...)
           standardGeneric("export.bed15"))

setMethod("export.bed15", "ANY",
          function(object, con, expNames, ...)
          {
            export.ucsc(object, con, "bed15", expNames = expNames, ...)
          })

setMethod("export.bed15", "UCSCData",
          function(object, con, expNames, ...)
          {
            ## ensure we do not override existing track line parameter
            if (missing(expNames) && is(object@trackLine, "Bed15TrackLine"))
              expNames <- object@trackLine@expNames
            export.ucsc(object, con, "bed15", expNames = expNames, ...)
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
        line@expScale <- as.numeric(vals["expScale"])
        line@expStep <- as.numeric(vals["expStep"])
        line@expNames <- strsplit(vals["expNames"], ",", fixed=TRUE)[[1]]
        line
      })

setGeneric("import.bedGraph",
           function(con, genome = "hg18", asRangedData = TRUE, ...)
           standardGeneric("import.bedGraph"))

setMethod("import.bedGraph", "ANY",
          function(con, genome, asRangedData = TRUE)
          {
            import.ucsc(con, "bedGraph", TRUE, genome = genome,
                        asRangedData = asRangedData)
          })

setGeneric("import.bedGraphLines",
           function(con, genome = "hg18", asRangedData = TRUE, ...)
           standardGeneric("import.bedGraphLines"))

setMethod("import.bedGraphLines", "ANY",
          function(con, genome, asRangedData = TRUE) {
            import.bed(con, variant = "bedGraph", genome = genome,
                       asRangedData = asRangedData)
          })

setGeneric("export.bedGraph",
           function(object, con, ...)
           standardGeneric("export.bedGraph"))

setMethod("export.bedGraph", "ANY",
          function(object, con, ...)
          {
            export.ucsc(object, con, "bedGraph", ...)
          })

setGeneric("export.bedGraphLines",
           function(object, con, ...)
           standardGeneric("export.bedGraphLines"))

setMethod("export.bedGraphLines", "ANY",
          function(object, con, ...)
          {
            export.bed(object, con, "bedGraph", ...)
          })
