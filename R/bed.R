# Import/export of Browser Extended Display (BED) data

setGeneric("export.bed",
           function(object, con, variant = c("base", "wig", "bed15"),
                    color = NULL, ...)
           standardGeneric("export.bed"))


setMethod("export.bed", "ANY",
          function(object, con, variant = c("base", "wig", "bed15"), color)
          {
            cl <- class(object)
            object <- try(as(object, "RangedData"), silent = TRUE)
            if (class(object) == "try-error")
              stop("cannot export object of class '", cl, "'")
            export.bed(object, con=con, variant=variant, color=color)
          })

setMethod("export.bed", c("RangedData", "characterORconnection"),
          function(object, con, variant = c("base", "wig", "bed15"), color)
          {
            variant <- match.arg(variant)
            name <- strand <- thickStart <- thickEnd <- color <- NULL
            blockCount <- blockSizes <- blockStarts <- NULL
            df <- data.frame(chrom(object), start(object)-1, end(object))
            score <- score(object)
            if (!is.null(score)) {
              if (!is.numeric(score) || any(is.na(score)))
                stop("Scores must be non-NA numeric values")
              if (variant != "wig" && any(score < 0 | score > 1000))
                stop("BED requires scores to fall within [0, 1000]")
            }
            if (variant == "wig") {
              if (is.null(score)) ## wig requires score
                score <- 0
              df$score <- score
            } else {
              blockSizes <- object$blockSizes
              blockStarts <- object$blockStarts
              if (variant == "bed15" && is.null(blockSizes))
                blockStarts <- blockSizes <- "" # bed15 must have all 15 cols
              if (!is.null(blockSizes))
                blockCount <- length(strsplit(blockSizes, ",")[[1]])
              if (is.null(color))
                color <- object$color
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
              df$color <- color
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
                        row.names = FALSE, quote = FALSE, na = ".")
          })

setMethod("export.bed", "UCSCData",
          function(object, con, variant = c("base", "wig", "bed15"), color,
                   trackLine = TRUE, ...)
          {
            variant <- match.arg(variant)
            if (variant == "base" && trackLine) {
              export.ucsc(object, con, "bed", variant, color, ...)
            } else {
              callNextMethod(object, con, variant, color)
            }
          })
          
setGeneric("import.bed",
           function(con, variant = c("base", "wig", "bed15"),
                    trackLine = TRUE, genome = "hg18", ...)
           standardGeneric("import.bed"))

setMethod("import.bed", "connection",
          function(con, variant = c("base", "wig", "bed15"), trackLine, genome)
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
                                   trackLine = FALSE, genome = genome))
            }
            if (variant == "wig") {
              bedClasses <- c("factor", "integer", "integer", "numeric")
              bedNames <- c("chrom", "start", "end", "score")
            } else {
              bedNames <- c("chrom", "start", "end", "name",
                            "score", "strand", "thickStart",
                            "thickEnd", "color", "blockCount",
                            "blockSizes", "blockStarts")
              bedClasses <- NA
            }
            if (variant == "bed15")
              bedNames <- c(bedNames, "expCount", "expIds", "expScores")
            ## FIXME: could read a single line to get ncols up-front,
            ## and thus specify all col classes
            ## FIXME: reading in 'as.is' to save memory,
            bed <- DataFrame(read.table(con, colClasses = bedClasses,
                                        as.is = TRUE))
            colnames(bed) <- bedNames[seq_len(ncol(bed))]
            if (variant != "wig") { ## don't know how many columns, coerce here
              bed$start <- as.integer(bed$start)
              bed$end <- as.integer(bed$end)
            } ## BED is 0-start, so add 1 to start
            color <- bed$color
            if (is.character(color)) { # could be NULL
              spec <- color != "0"
              cols <- unlist(strsplit(color[spec], ",", fixed=TRUE),
                             use.names=FALSE)
              cols <- matrix(as.integer(cols), 3)
              color <- rep(NA, nrow(bed))
              color[spec] <- rgb(cols[1,], cols[2,], cols[3,], max = 255)
              bed$color <- color              
            }
            GenomicData(IRanges(bed$start + 1, bed$end),
                        bed[,tail(colnames(bed), -3),drop=FALSE],
                        chrom = bed$chrom, genome = genome)
          })

setGeneric("blocks", function(x, ...) standardGeneric("blocks"))

setMethod("blocks", "RangedData",
          function(x)
          {
            if (is.null(x$blockStarts) || is.null(x$blockSizes))
              stop("'x' must have 'blockStarts' and 'blockSizes' columns")
            starts <- unlist(strsplit(as.character(x$blockStarts), ","))
            starts <- as.integer(starts) + 1
            sizes <- unlist(strsplit(as.character(x$blockSizes), ","))
            split(IRanges(starts, width = sizes), x$name)
          })

setGeneric("import.bed15Lines",
           function(con, trackLine, genome = "hg18", ...)
           standardGeneric("import.bed15Lines"))

setMethod("import.bed15Lines", "ANY",
          function(con, trackLine, genome)
          {
            bed <- import.bed(con, "bed15", genome = genome)
            if (!nrow(bed))
              return(bed)
            ids <- strsplit(bed$expIds[1], ",", fixed=TRUE)[[1]]
            expNames <- trackLine@expNames[as.integer(ids) + 1]
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
          })

setGeneric("import.bed15",
           function(con, genome = "hg18", ...) standardGeneric("import.bed15"))

setMethod("import.bed15", "ANY",
          function(con, genome)
          {
            import.ucsc(con, "bed15", TRUE, genome = genome)
          })

setGeneric("export.bed15Lines",
           function(object, con, trackLine, ...)
           standardGeneric("export.bed15Lines"))


setMethod("export.bed15Lines", "ANY",
          function(object, con, trackLine, ...)
          {
            cl <- class(object)
            object <- try(as(object, "RangedData"), silent = TRUE)
            if (class(object) == "try-error")
              stop("cannot export object of class '", cl, "'")
            export.bed15Lines(object, con=con, trackLine=trackLine, ...)
          })

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
            export.bed(object, con, "bed15")
          })

setGeneric("export.bed15",
           function(object, con, expNames = NULL, ...)
           standardGeneric("export.bed15"))

setMethod("export.bed15", "ANY",
          function(object, con, expNames, ...)
          {
            cl <- class(object)
            object <- try(as(object, "RangedData"), silent = TRUE)
            if (class(object) == "try-error")
              stop("cannot export object of class '", cl, "'")
            export.bed15(object, con=con, expNames=expNames, ...)
          })

setMethod("export.bed15", "RangedData",
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
