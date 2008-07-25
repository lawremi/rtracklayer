# import/export of wiggle (WIG) format

setGeneric("import.wig", function(con, ...) standardGeneric("import.wig"))
setMethod("import.wig", "ANY",
          function(con)
          {
            import.ucsc(con, "wig")
          })

setGeneric("export.wig",
           function(object, con,
                    dataFormat = c("auto", "bed", "variableStep", "fixedStep"),
                    ...)
           standardGeneric("export.wig"))
setMethod("export.wig", "trackSet",
          function(object, con,
                   dataFormat = c("auto", "bed", "variableStep", "fixedStep"),
                   ...)
          {
            export.ucsc(object, con, "wig", dataFormat = match.arg(dataFormat),
                        ...)
          })

setGeneric("export.wigLines",
           function(object, con,
                    dataFormat = c("auto", "bed", "variableStep", "fixedStep"),
                    ...)
           standardGeneric("export.wigLines"))
setMethod("export.wigLines", "trackSet",
          function(object, con,
                   dataFormat = c("auto", "bed", "variableStep", "fixedStep"))
          {
            vals <- sampleNames(object)[1] # can only use a single condition
            object <- object[!is.na(dataVals(object)[,1])]
            object <- object[order(start(object))]
            object <- object[order(chrom(object))]
            df <- trackData(object)
            ## df <- df[!is.na(df[[vals]]),] # no NAs
            ## df <- df[order(df$start),]
            ## attempt to use most efficient format if not specified
### FIXME: If we need bed, use bedGraph format
            byChrom <- function(chromData, formatOnly) {
              starts <- chromData$start
              ends <- chromData$end
              if (!all(tail(starts, -1) - head(ends, -1) > 0))
                stop("Features cannot overlap. ",
                     "Note that WIG does not distinguish between strands - ",
                     "try exporting two tracks, one for each strand.")
              spans <- ends - starts + 1
              if (length(starts) == 1)
                steps <- 0
              else steps <- diff(starts)
              fixedSpan <- all(spans[1] == spans)
              fixedStep <- all(steps[1] == steps)
              dataFormat <- "bed"
              if (fixedSpan) {
                dataFormat <- "variableStep"
                if (fixedStep)
                  dataFormat <- "fixedStep"
              }
              if (formatOnly)
                return(dataFormat)
              cat(dataFormat, file = con)
              cat(" chrom=", as.character(chromData$chrom)[1],
                  file = con, sep = "")
              data <- chromData[[vals]]
              if (dataFormat == "variableStep")
                data <- cbind(starts, data)
              else {
                if (!fixedStep)
                  stop("Format 'fixedStep' invalid: step not uniform")
                cat(" start=", starts[1], file = con, sep = "")
                cat(" step=", steps[1], file = con, sep = "")
              }
              if (fixedSpan)
                cat(" span=", spans[1], file = con, sep = "")
              cat("\n", file = con, sep = "")
              write.table(data, con, sep = "\t", col.names = FALSE,
                          row.names = FALSE, quote = FALSE)
            }
            dataFormat <- match.arg(dataFormat)
            if (dataFormat == "auto")
              dataFormat <- by(df, df$chrom, byChrom, TRUE)
            if (any(dataFormat == "bed")) # one BED, all BED
              export.bed(object, con, wig = TRUE)
            else by(df, df$chrom, byChrom, FALSE) # else, mix variable/fixed
          })

setGeneric("import.wig",
           function(con, genome = "hg18", ...) standardGeneric("import.wig"))
setMethod("import.wig", "ANY",
          function(con, genome)
          {
            import.ucsc(con, "wig", TRUE, genome = genome)
          })

setGeneric("import.wigLines",
           function(con, genome = "hg18", ...) standardGeneric("import.wigLines"))
setMethod("import.wigLines", "ANY",
          function(con, genome)
          {
            lines <- readLines(con, warn = FALSE)
            formatInds <- grep("^variableStep|^fixedStep", lines)
            formatLines <- lines[formatInds]
            starts <- formatInds+1
            ends <- c(tail(formatInds,-1)-1, length(lines))
            if (length(formatLines)) {
              parseData <- function(i) {
                # parse the data values
                con <- file()
                writeLines(lines[starts[i]:ends[i]], con)
                data <- read.table(con)
                close(con)
                # parse format line
                format <- gsub("^([^ ]*) .*", "\\1", formatLines[i])
                formatVals <- ucscParsePairs(formatLines[i])
                if (format == "variableStep") {
                  start <- data[,1]
                  score <- data[,2]
                } else {
                  start <- seq(as.numeric(formatVals["start"]),
                               by = as.numeric(formatVals["step"]),
                               length.out = nrow(data))
                  score <- data[,1]
                }
                span <- formatVals["span"]
                if (is.na(span))
                  span <- 1
                end <- start + as.numeric(span) - 1
                data.frame(chrom = formatVals[["chrom"]],
                           start = start, end = end,
                           score = score)
              }
              resultList <- lapply(seq_along(formatInds), parseData)
              resultMat <- do.call("rbind", resultList)
              featureData <- resultMat[,c("chrom", "start", "end")]
              new("trackSet", featureData = trackFeatureData(featureData),
                  dataVals = resultMat[,"score",drop=FALSE], genome = genome)
            } else import(text = lines, format = "bed", wig = TRUE)
        })

setClass("wigTrackLine",
         representation(altColor = "integer", autoScale = "logical",
                        gridDefault = "logical", maxHeightPixels = "numeric",
                        graphType = "character", viewLimits = "numeric",
                        yLineMark = "numeric", yLineOnOff = "logical",
                        windowingFunction = "character",
                        smoothingWindow = "numeric"),
         contains = "ucscTrackLine")

setAs("wigTrackLine", "character",
      function(from)
      {
        str <- as(as(from, "ucscTrackLine"), "character")
        str <- paste(str, "type=wiggle_0")
        color <- from@altColor
        if (length(color))
          str <- paste(str, " altColor=", paste(color, collapse=","), sep="")
        autoScale <- from@autoScale
        if (length(autoScale) && !autoScale)
          str <- paste(str, "autoScale=off")
        gridDefault <- from@gridDefault
        if (length(gridDefault) && gridDefault)
          str <- paste(str, "gridDefault=On")
        maxHeightPixels <- from@maxHeightPixels
        if (length(maxHeightPixels) && maxHeightPixels)
          str <- paste(str, " maxHeightPixels=",
                       paste(maxHeightPixels, collapse=":"), sep = "")
        graphType <- from@graphType
        if (length(graphType))
          str <- paste(str, " graphType=", graphType, sep = "")
        viewLimits <- from@viewLimits
        if (length(viewLimits))
          str <- paste(str, " viewLimits=", paste(viewLimits, collapse = ":"),
                       sep = "")
        yLineMark <- from@yLineMark
        if (length(yLineMark))
          str <- paste(str, " yLineMark=", yLineMark, sep = "")
        yLineOnOff <- from@yLineOnOff
        if (length(yLineOnOff) && yLineOnOff)
          str <- paste(str, "yLineOnOff=On")
        windowingFunction <- from@windowingFunction
        if (length(windowingFunction))
          str <- paste(str, " windowingFunction=", windowingFunction, sep = "")
        smoothingWindow <- from@smoothingWindow
        if (length(smoothingWindow))
          str <- paste(str, " smoothingWindow=", smoothingWindow, sep = "")
        str
      })

setAs("character", "wigTrackLine",
      function(from)
      {
        line <- new("wigTrackLine", as(from, "ucscTrackLine"))
        vals <- ucscParsePairs(from)
        if (!is.na(vals["altColor"]))
          line@altColor <- as.integer(strsplit(vals["altColor"], ",")[[1]])
        if (!is.na(vals["autoScale"]))
          line@autoScale <- vals["autoScale"] == "On"
        if (!is.na(vals["gridDefault"]))
          line@gridDefault <- vals["gridDefault"] == "On"
        if (!is.na(vals["maxHeightPixels"]))
          line@maxHeightPixels <-
            as.numeric(strsplit(vals["maxHeightPixels"], ":")[[1]])
        if (!is.na(vals["graphType"]))
          line@graphType <- vals["graphType"]
        if (!is.na(vals["viewLimits"]))
          line@viewLimits <-
            as.numeric(strsplit(vals["viewLimits"], ":")[[1]])
        if (!is.na(vals["yLineMark"]))
          line@yLineMark <- as.numeric(vals["yLineMark"])
        if (!is.na(vals["yLineOnOff"]))
          line@yLineOnOff <- vals["yLineOnOff"] == "On"
        if (!is.na(vals["windowingFunction"]))
          line@windowingFunction <- vals["windowingFunction"]
        if (!is.na(vals["smoothingWindow"]))
          line@smoothingWindow <- as.numeric(vals["smoothingWindow"])
        line
      })

setAs("basicTrackLine", "wigTrackLine",
      function(from) new("wigTrackLine", from))

setAs("wigTrackLine", "basicTrackLine",
      function(from) new("basicTrackLine", from))
