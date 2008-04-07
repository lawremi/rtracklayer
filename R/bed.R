# Import/export of Browser Extended Display (BED) data

setGeneric("export.bed",
           function(object, con, wig = FALSE, ...)
           standardGeneric("export.bed"))

setMethod("export.bed", "trackSet",
          function(object, con, wig)
          {
            df <- trackData(object)
            bed <- cbind(as.character(df$featChrom), df$featStart, df$featEnd)
            if (!wig)
              bed <- cbind(bed, featureNames(object))
            bed <- cbind(bed, df[[sampleNames(object)[1]]])
            if (!wig) {
              blockCount <- NULL
              if (!is.null(df$blockSizes))
                blockCount <- length(strsplit(df$blockSizes, ",")[[1]])
              color <- df$color
              if (is.null(color) && !is.null(blockCount))
                color <- "black"
              if (!is.null(color))
                color <- col2rgb(color)[,1]
              thickStart <- df$thickStart
              thickEnd <- df$thickEnd
              if (is.null(thickStart) && !is.null(color)) {
                thickStart <- df$featStart
                thickEnd <- df$featEnd
              }
              bed <- cbind(bed, df$featStrand, thickStart, thickEnd, color,
                           blockCount, df$blockSizes, df$blockStarts)
            }
            write.table(bed, con, sep = "\t", col.names = FALSE,
                        row.names = FALSE, quote = FALSE)
          })

setMethod("export.bed", "ucscTrackSet",
          function(object, con, wig, trackLine = !wig)
          {
            if (!wig && trackLine)
              export.ucsc(object, con, "bed", wig)
            else export.wig(as(object, "trackSet"), con, wig)
          })
          
setGeneric("import.bed",
           function(con, wig = FALSE, trackLine = !wig, ...)
           standardGeneric("import.bed"))

setMethod("import.bed", "ANY",
          function(con, wig, trackLine)
          {
            if (!wig && trackLine) {
              # check for a track line
              lines <- readLines(con, warn = FALSE)
              if (length(grep("^track", lines)) > 0)
                trackSet <- import(text = lines, format = "ucsc",
                                   subformat = "bed", drop = TRUE,
                                   trackLine = FALSE)
              else trackLine <- FALSE
            }
            if (wig || !trackLine) {
              bed <- read.table(con)
              bedNames <- c("featChrom", "featStart", "featEnd", "name",
                            "score", "featStrand", "thickStart",
                            "thickEnd", "color", "blockCount", "blockSizes",
                            "blockStarts")
              colnames(bed) <- bedNames[seq_len(ncol(bed))]
              featureData <- bed[,!(colnames(bed) %in% c("name", "score"))]
              if (!wig)
                rownames(featureData) <- bed$name
              score <- NA
              if (wig)
                score <- bed$name
              else if (!is.null(bed$score))
                score <- bed$score
              trackSet <- new("trackSet", featureData = featureData,
                              dataVals = score)
            }
            trackSet
          })
