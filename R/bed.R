# Import/export of Browser Extended Display (BED) data

setGeneric("export.bed",
           function(object, con, wig = FALSE, color = NULL, ...)
           standardGeneric("export.bed"))

setMethod("export.bed", "trackSet",
          function(object, con, wig, color)
          {
            df <- trackData(object) # BED start positions are 0-based
            bed <- cbind(as.character(df$chrom), as.numeric(df$start) - 1,
                         df$end)
            if (!wig) {
              name <- as.character(df$name)
              if (!length(name))
                name <- featureNames(object)
              bed <- cbind(bed, name)
            }
            score <- df[[sampleNames(object)[1]]]
            if (is.null(score))
              score <- 0
            bed <- cbind(bed, score)
            if (!wig) {
              blockCount <- NULL
              if (!is.null(df$blockSizes))
                blockCount <- length(strsplit(df$blockSizes, ",")[[1]])
              if (is.null(color))
                color <- df$color
              if (is.null(color) && !is.null(blockCount))
                color <- "black"
              if (!is.null(color))
                color <- paste(as.vector(col2rgb(color)), collapse=",")
              thickStart <- df$thickStart
              thickEnd <- df$thickEnd
              if (is.null(thickStart) && !is.null(color)) {
                thickStart <- df$start
                thickEnd <- df$end
              }
              bed <- cbind(bed, as.character(df$strand), thickStart, thickEnd,
                           color, blockCount, df$blockSizes, df$blockStarts)
            }
            write.table(bed, con, sep = "\t", col.names = FALSE,
                        row.names = FALSE, quote = FALSE, na = ".")
          })

setMethod("export.bed", "ucscTrackSet",
          function(object, con, wig, trackLine = !wig, color, ...)
          {
            if (!wig && trackLine) {
              export.ucsc(object, con, "bed", wig, color, ...)
            } else export.bed(as(object, "trackSet"), con, wig, color)
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
                trackSet <- import.ucsc(con, subformat = "bed", drop = TRUE,
                                        trackLine = FALSE, genome = genome)
              else trackLine <- FALSE
            }
            if (wig || !trackLine) {
              bed <- read.table(con)
              bedNames <- c("chrom", "start", "end", "name",
                            "score", "strand", "thickStart",
                            "thickEnd", "color", "blockCount", "blockSizes",
                            "blockStarts")
              colnames(bed) <- bedNames[seq_len(ncol(bed))]
              bed$start <- bed$start + 1 # BED has 0-based start positions
              featureData <- bed[,!(colnames(bed) == "score")]
              if (!wig)
                rownames(featureData) <- make.names(bed$name, TRUE)
              if (wig) {
                score <- bed$name
                featureData$name <- NULL
              } else if (!is.null(bed$score))
                score <- bed$score
              else score <- rep(NA, nrow(featureData))
              trackSet <- new("trackSet",
                              featureData = trackFeatureData(featureData),
                              dataVals = cbind(score = score), genome = genome)
            }
            trackSet
          })
