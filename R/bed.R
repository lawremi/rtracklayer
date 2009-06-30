# Import/export of Browser Extended Display (BED) data

setGeneric("export.bed",
           function(object, con, wig = FALSE, color = NULL, ...)
           standardGeneric("export.bed"))

setMethod("export.bed", "RangedData",
          function(object, con, wig, color)
          {
            name <- strand <- thickStart <- thickEnd <- color <- NULL
            blockCount <- blockSizes <- blockStarts <- NULL
            df <- data.frame(chrom(object), start(object)-1, end(object))
            score <- score(object)
            if (!is.null(score)) {
              if (!is.numeric(score) || any(is.na(score)))
                stop("Scores must be non-NA numeric values")
              if (!wig && any(score < 0 | score > 1000))
                stop("BED requires scores to fall within [0, 1000]")
            }
            if (wig) {
              if (is.null(score)) ## wig requires score
                score <- 0
              df$score <- score
            } else {
              blockSizes <- object$blockSizes
              if (!is.null(blockSizes)) {
                blockCount <- length(strsplit(blockSizes, ",")[[1]])
                blockStarts <- object$blockStarts
              }
              if (is.null(color))
                color <- object$color
              if (is.null(color) && !is.null(blockCount))
                color <- "black" ## blocks require color
              if (!is.null(color))
                color <- paste(as.vector(col2rgb(color)), collapse=",")
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
            }
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            write.table(df, con, sep = "\t", col.names = FALSE,
                        row.names = FALSE, quote = FALSE, na = ".")
          })

setMethod("export.bed", "UCSCData",
          function(object, con, wig, color, trackLine = !wig, ...)
          {
            if (!wig && trackLine) {
              export.ucsc(object, con, "bed", wig, color, ...)
            } else callNextMethod(object, con, wig, color)
          })
          
setGeneric("import.bed",
           function(con, wig = FALSE, trackLine = !wig, genome = "hg18", ...)
           standardGeneric("import.bed"))

setMethod("import.bed", "ANY",
          function(con, wig, trackLine, genome)
          {
            if (!wig && trackLine) {
              ## check for a track line
              line <- "#"
              while(length(grep("^ *#", line))) # skip initial comments
                line <- readLines(con, 1, warn = FALSE)
              pushBack(line, con)
              if (length(grep("^track", line)) > 0)
                return(import.ucsc(con, subformat = "bed", drop = TRUE,
                                   trackLine = FALSE, genome = genome))
            }
            if (wig) {
              bedClasses <- c("factor", "integer", "integer", "numeric")
              bedNames <- c("chrom", "start", "end", "score")
            } else {
              bedNames <- c("chrom", "start", "end", "name",
                            "score", "strand", "thickStart",
                            "thickEnd", "color", "blockCount",
                            "blockSizes", "blockStarts")
              bedClasses <- NA
            }
            ## FIXME: could read a single line to get ncols up-front,
            ## and thus specify all col classes
            ## FIXME: reading in 'as.is' to save memory,
            ## XFactor would be useful here.
            bed <- DataFrame(read.table(con, colClasses = bedClasses,
                                        as.is = TRUE))
            colnames(bed) <- bedNames[seq_len(ncol(bed))]
            if (!wig) { ## don't know how many columns, coerce here
              bed$start <- as.integer(bed$start)
              bed$end <- as.integer(bed$end)
            } ## BED is 0-start, so add 1 to start
            track <- GenomicData(IRanges(bed$start + 1, bed$end),
                                 bed[,tail(colnames(bed), -3),drop=FALSE],
                                 chrom = bed$chrom, genome = genome)
            
            ##featureData <- bed[,!(colnames(bed) == "score")]
            ##if (!wig)
            ##  rownames(featureData) <- make.names(bed$name, TRUE)
            ##if (wig) {
            ##  score <- bed$name
            ##  featureData$name <- NULL
            ##} else if (!is.null(bed$score))
            ##score <- bed$score
            ##else score <- rep(NA, nrow(featureData))
            ##trackSet <- new("trackSet",
            ##                featureData = trackFeatureData(featureData),
            ##               dataVals = cbind(score = score), genome = genome)
            track
          })

setGeneric("blocks", function(x, ...) standardGeneric("blocks"))

setMethod("blocks", "RangedData",
          function(x)
          {
            if (is.null(x$blockStarts) || is.null(x$blockSizes))
              stop("'x' must have 'blockStarts' and 'blockSizes' columns")
            starts <- unlist(strsplit(as.character(x$blockStarts), ",")) + 1
            sizes <- unlist(strsplit(as.character(x$blockSizes), ","))
            split(IRanges(starts, width = sizes), x$name)
          })
