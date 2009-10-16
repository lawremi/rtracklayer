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
            export.ucsc(object, con, "wig", dataFormat = match.arg(dataFormat),
                        ...)
          })

setGeneric("export.wigLines",
           function(object, con,
                    dataFormat = c("auto", "variableStep", "fixedStep"),
                    ...)
           standardGeneric("export.wigLines"))

setMethod("export.wigLines", c("RangedData", "characterORconnection"),
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep"))
          {
            if (any(is.na(score(object))))
              stop("WIG cannot encode missing values")
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            doBlock <- function(chromData) {
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
              ## heuristic: split into blocks with fixed span
              fixedSpan <- all(spans[1] == spans)
              fixedStep <- all(steps[1] == steps)
              if (dataFormat == "auto") {
                dataFormat <- "variableStep"
                if (fixedStep)
                  dataFormat <- "fixedStep"
              }
              if (!fixedSpan) ## split into blocks to make spans uniform
                return(lapply(split(chromData, spans), doBlock))
              cat(dataFormat, file = con)
              cat(" chrom=", as.character(space(chromData))[1],
                  file = con, sep = "")
              data <- score(chromData)
              if (dataFormat == "variableStep")
                data <- cbind(starts, data)
              else {
                if (!fixedStep)
                  stop("Step not uniform: consider variableStep")
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
            lapply(object, doBlock)
          })

setGeneric("import.wig",
           function(con, genome = "hg18", ...) standardGeneric("import.wig"))
setMethod("import.wig", "ANY",
          function(con, genome)
          {
            import.ucsc(con, "wig", TRUE, genome = genome)
          })

setGeneric("import.wigLines",
           function(con, genome = "hg18", ...)
           standardGeneric("import.wigLines"))

setMethod("import.wigLines", "characterORconnection",
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
                            chrom = formatVals[["chrom"]])
              }
              resultList <- lapply(seq_along(formatInds), parseData)
              gd <- do.call(rbind, resultList)
              genome(gd) <- genome
              gd
            } else import(text = lines, format = "bed", variant = "bedGraph",
                          genome = genome)
        })
