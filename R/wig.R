# import/export of wiggle (WIG) format

setGeneric("export.wig",
           function(object, con,
                    dataFormat = c("auto", "variableStep", "fixedStep"),
                    ...)
           standardGeneric("export.wig"))
setMethod("export.wig", "ANY",
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep"),
                   ...)
          {
            export.ucsc(object, con, "wig", 
                        dataFormat = match.arg(dataFormat), ...)
          })

setGeneric("export.wigLines",
           function(object, con,
                    dataFormat = c("auto", "variableStep", "fixedStep"),
                    ...)
           standardGeneric("export.wigLines"))

.wigWriter <- function(chromData, con, dataFormat) {
  cat(dataFormat, file = con, append = TRUE)
  cat(" chrom=", as.character(space(chromData))[1],
      file = con, sep = "", append = TRUE)
  data <- score(chromData)
  starts <- start(chromData)
  if (dataFormat == "variableStep")
    data <- cbind(starts, data)
  else {
    cat(" start=", starts[1], file = con, sep = "", append = TRUE)
    step <- if (length(starts) == 1) 0 else starts[2] - starts[1]
    cat(" step=", step, file = con, sep = "", append = TRUE)
  }  
  cat(" span=", width(chromData)[1], file = con, sep = "", append = TRUE)
  cat("\n", file = con, sep = "", append = TRUE)
  write.table(data, con, sep = "\t", col.names = FALSE,
              row.names = FALSE, quote = FALSE, append = TRUE)
}

setMethod("export.wigLines", c("RangedData", "characterORconnection"),
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep"),
                   writer = .wigWriter, append = FALSE)
          {
            if (!is.numeric(score(object)) || any(is.na(score(object))))
              stop("The score must be numeric, without any NA's")
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            if (!append)
              cat("", file = con)
            dataFormat <- match.arg(dataFormat)            
            doBlock <- function(chromData) {
              if (!nrow(chromData))
                return()
              if (is.unsorted(start(chromData)))
                chromData <- chromData[order(start(chromData)),]
              starts <- start(chromData)
              ends <- end(chromData)
              if (!all(tail(starts, -1) - head(ends, -1) > 0))
                stop("Features cannot overlap. ",
                     "Note that WIG does not distinguish between strands - ",
                     "try exporting two tracks, one for each strand.")
              spans <- ends - starts + 1
              if (length(starts) == 1)
                steps <- 0
              else steps <- diff(starts)
              fixedSpan <- all(spans[1] == spans)
              if (!fixedSpan)
                stop("The span must be uniform for Wiggle export. ",
                     "Consider bedGraph or bigWig as alternatives.")
              fixedStep <- all(steps[1] == steps)
              if (dataFormat == "auto") {
                dataFormat <- "variableStep"
                if (fixedStep)
                  dataFormat <- "fixedStep"
              }
              if (dataFormat != "variableStep" && !fixedStep)
                stop("Step not uniform: consider variableStep")
              writer(chromData, con, dataFormat)
            }
            dataFormat <- match.arg(dataFormat)
            invisible(lapply(object, doBlock))
          })

setGeneric("import.wig",
           function(con, genome = NULL, asRangedData = TRUE, ...)
           standardGeneric("import.wig"))
setMethod("import.wig", "ANY",
          function(con, genome, asRangedData = TRUE)
          {
            import.ucsc(con, "wig", TRUE, genome = genome,
                        asRangedData = asRangedData)
          })

setGeneric("import.wigLines",
           function(con, genome = NULL, asRangedData = TRUE, ...)
           standardGeneric("import.wigLines"))

setMethod("import.wigLines", "characterORconnection",
          function(con, genome, asRangedData = TRUE)
          {
            if (!isTRUEorFALSE(asRangedData))
              stop("'asRangedData' must be TRUE or FALSE")
            lines <- readLines(con, warn = FALSE)
            formatInds <- grep("^variableStep|^fixedStep", lines)
            formatLines <- lines[formatInds]
            starts <- formatInds + 1L
            ends <- c(tail(formatInds, -1) - 1L, length(lines))
            if (length(formatLines)) {
              parseData <- function(i) {
                # parse the data values
                con <- file()
                writeLines(window(lines, starts[i], ends[i]), con)
                data <- read.table(con)
                close(con)
                # parse format line
                format <- gsub("^([^ ]*) .*", "\\1", formatLines[i])
                formatVals <- ucscParsePairs(formatLines[i])
                if (format == "variableStep") {
                  start <- data[,1]
                  score <- data[,2]
                } else {
                  start <- seq(as.integer(formatVals["start"]),
                               by = as.integer(formatVals["step"]),
                               length.out = nrow(data))
                  score <- data[,1]
                }
                span <- formatVals["span"]
                if (is.na(span))
                  span <- 1
                end <- start + as.integer(span) - 1
                GenomicData(IRanges(start, end), score = score,
                            chrom = formatVals[["chrom"]],
                            asRangedData = asRangedData)
              }
              resultList <- lapply(seq_along(formatInds), parseData)
              if (asRangedData)
                gd <- do.call(rbind, resultList)
              else
                gd <- do.call(c, resultList)
              genome(gd) <- genome
              gd
            } else {
              import(text = lines, format = "bed", variant = "bedGraph",
                     genome = genome, asRangedData = asRangedData)
            }
        })
