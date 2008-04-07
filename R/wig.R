# import/export of wiggle (WIG) format

setGeneric("import.wig", function(con, ...) standardGeneric("import.wig"))
setMethod("import.wig", "ANY",
          function(con)
          {
            import.ucsc(con, "wig")
          })

setGeneric("export.wig",
           function(object, con,
                    dataFormat = c("bed", "variableStep", "fixedStep"), ...)
           standardGeneric("export.wig"))
setMethod("export.wig", "trackSet",
          function(object, con, dataFormat)
          {
            export.ucsc(object, con, "wig", dataFormat)
          })

setGeneric("export.wigLines",
           function(object, con,
                    dataFormat = c("bed", "variableStep", "fixedStep"), ...)
           standardGeneric("export.wigLines"))
setMethod("export.wigLines", "trackSet",
          function(object, con,
                   dataFormat = c("bed", "variableStep", "fixedStep"))
          {
            formatMissing <- length(dataFormat) > 1
            formatMatched <- match.arg(dataFormat)
            vals <- sampleNames(object)[1] # can only use a single condition
            df <- trackData(object)
            by(df, df$featChrom,
               function(chromData)
               {
                 chromData <- chromData[!is.na(chromData[[vals]]),] # no NAs
                 chromData <- chromData[order(chromData$featStart),]
                 starts <- chromData$featStart
                 ends <- chromData$featEnd
                 spans <- ends - starts
                 steps <- diff(starts)
                 fixedSpan <- all(spans[1] == spans)
                 fixedStep <- all(steps[1] == steps)
                 # attempt to use most efficient format if not specified
                 if (formatMissing && fixedSpan) { # fixed span
                   dataFormat <- "variableStep"
                   if (fixedStep) # and a fixed step
                     dataFormat <- "fixedStep"
                 } else dataFormat <- formatMatched
                 if (dataFormat == "bed")
                   export.bed(object, con, wig = TRUE)
                 else {
                   cat(dataFormat, file = con)
                   cat(" chrom=", as.character(chromData$featChrom)[1],
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
               })
          })

setGeneric("import.wig",
           function(con, ...) standardGeneric("import.wig"))
setMethod("import.wig", "ANY",
          function(con)
          {
            import.ucsc(con, "wig", TRUE)
          })

setGeneric("import.wigLines",
           function(con, ...) standardGeneric("import.wigLines"))
setMethod("import.wigLines", "ANY",
          function(con)
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
                  featStart <- data[,1]
                  score <- data[,2]
                } else {
                  featStart <- seq(as.numeric(formatVals["start"]),
                                  by = as.numeric(formatVals["step"]),
                                  length.out = nrow(data))
                  score <- data[,1]
                }
                featEnd <- c(tail(featStart,-1)-1, tail(featStart, 1))
                if (!is.na(formatVals["span"]))
                  featEnd <- featStart + as.numeric(formatVals["span"])
                data.frame(featChrom = formatVals[["chrom"]],
                           featStart = featStart, featEnd = featEnd,
                           score = score)
              }
              resultList <- lapply(seq_along(formatInds), parseData)
              resultMat <- do.call("rbind", resultList)
              featureData <- resultMat[,c("featChrom", "featStart", "featEnd")]
              new("trackSet", featureData = featureData,
                  dataVals = resultMat[,"score"])
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
        if (length(autoScale) && autoScale)
          str <- paste(str, "autoScale=On")
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
