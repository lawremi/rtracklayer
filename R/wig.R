### =========================================================================
### WIG (Wiggle) support (fixedStep and variableStep, legacy bedGraph)
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("WIGFile", contains = "RTLFile")

WIGFile <- function(resource) {
  new("WIGFile", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.wig",
           function(object, con,
                    dataFormat = c("auto", "variableStep", "fixedStep"),
                    ...)
           standardGeneric("export.wig"),
           signature = c("object", "con"))

setMethod("export.wig", "ANY",
          function(object, con,
                   dataFormat = c("auto", "variableStep", "fixedStep"),
                   ...)
          {
            export(object, con, "wig", dataFormat = match.arg(dataFormat), ...)
          })

setMethod("export", c("ANY", "WIGFile"),
          function(object, con, ...)
          {
            cl <- class(object)
            track <- try(as(object, "RangedData"), silent = TRUE)
            if (class(track) == "try-error") {
              track <- try(as(object, "RangedDataList"), silent = TRUE)
              if (class(track) == "try-error")
                stop("cannot export object of class '", cl, "'")
            }
            export(track, con, ...)
          })

.wigWriter <- function(chromData, con, dataFormat, append) {
  con <- connection(con, if (append) "a" else "w")
  on.exit(release(con))
  cat(dataFormat, file = con)
  cat(" chrom=", as.character(space(chromData))[1], file = con, sep = "")
  data <- score(chromData)
  starts <- start(chromData)
  if (dataFormat == "variableStep")
    data <- cbind(starts, data)
  else {
    cat(" start=", starts[1], file = con, sep = "")
    step <- if (length(starts) == 1) 0 else starts[2] - starts[1]
    cat(" step=", step, file = con, sep = "")
  }  
  cat(" span=", width(chromData)[1], file = con, sep = "")
  cat("\n", file = con, sep = "")
  write.table(data, con, sep = "\t", col.names = FALSE,
              row.names = FALSE, quote = FALSE)
}

setMethod("export", c("RangedData", "WIGFile"),
          function(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep"),
                   writer = .wigWriter, append = FALSE, trackLine = TRUE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (trackLine) {
              return(export.ucsc(object, con,
                                 dataFormat = match.arg(dataFormat),
                                 trackLine = FALSE, append = append))
            }
            if (!is.numeric(score(object)) || any(is.na(score(object))))
              stop("The score must be numeric, without any NA's")
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
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
              ans <- writer(chromData, con, dataFormat, append)
              append <<- TRUE
              ans
            }
            dataFormat <- match.arg(dataFormat)
            invisible(lapply(object, doBlock))
          })

setMethod("export", c("RangedDataList", "WIGFile"),
          .export_RangedDataList_RTLFile)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.wig",
           function(con, genome = NA, asRangedData = TRUE, which = NULL, ...)
           standardGeneric("import.wig"),
           signature = "con")

setMethod("import.wig", "ANY",
          function(con, genome, asRangedData = TRUE, ...)
          {
            import(con, "wig", genome = genome, asRangedData = asRangedData,
                   which = which, ...)
          })

setMethod("import", "WIGFile",
          function(con, format, text, genome = NA, asRangedData = TRUE,
                   trackLine = TRUE, which = NULL)
          {
            file <- con
            con <- connection(con, "r")
            ## check for a track line
            line <- scanTrackLine(con)
            if (!is.null(line) && trackLine) {
              pushBack(line, con)
              return(import.ucsc(initialize(file, resource = con), drop = TRUE,
                                 trackLine = FALSE, genome = genome,
                                 asRangedData = asRangedData))
            }
            if (!isTRUEorFALSE(asRangedData))
              stop("'asRangedData' must be TRUE or FALSE")
            lines <- readLines(con, warn = FALSE)
            formatInds <- grep("chrom=", lines, fixed=TRUE)
            formatLines <- lines[formatInds]
            starts <- formatInds + 1L
            ends <- c(tail(formatInds, -1) - 1L, length(lines))
            format <- gsub("^([^ ]*) .*", "\\1", formatLines)
            parsedFormat <- lapply(formatLines, ucscParsePairs)
            chrom <- sapply(parsedFormat, `[[`, "chrom")
            seqlevels <- unique(chrom)
            if (length(formatLines)) {
              parseData <- function(i) {
                ## parse format line
                formatVals <- parsedFormat[[i]]
                # parse the data values
                block_lines <- window(lines, starts[i], ends[i])
                if (!length(block_lines)) {
                  if (asRangedData) {
                    rl <- RangesList(IRanges())
                    names(rl) <- formatVals[["chrom"]]
                    return(RangedData(rl, score = numeric()))
                  } else {
                    gr <- GRanges(score = numeric())
                    seqlevels(gr) <- seqlevels
                    return(gr)
                  }
                }
                block_lines <- grep("^#", block_lines, invert = TRUE,
                                    value = TRUE)
                if (format[i] == "variableStep") {
                  ## assume the same white space for every row
                  sep <- sub(".*?([[:space:]]+).*", "\\1", block_lines[1])
                  split_lines <- strsplit(block_lines, sep, fixed=TRUE)
                  mat <- matrix(as.numeric(unlist(split_lines)), nrow = 2)
                  start <- mat[1,]
                  score <- mat[2,]
                } else {
                  start <- seq(as.integer(formatVals["start"]),
                               by = as.integer(formatVals["step"]),
                               length.out = length(block_lines))
                  score <- as.numeric(block_lines)
                }
                span <- formatVals["span"]
                if (is.na(span))
                  span <- 1
                end <- start + as.integer(span) - 1
                gd <- GenomicData(IRanges(start, end), score = score,
                                  chrom = formatVals[["chrom"]],
                                  asRangedData = asRangedData,
                                  which = which)
                seqlevels(gd) <- seqlevels
                gd
              }
              parseInds <- seq_along(formatInds)
              if (!is.null(which)) {
                message("For efficiency, consider converting this WIG file",
                        " to a BigWig file; see ?wigToBigWig")
                parseInds <- parseInds[chrom %in% seqnames(which)]
              }
              resultList <- lapply(parseInds, parseData)
              if (asRangedData) {
                rl <- do.call(c, lapply(resultList, ranges))
                gd <- RangedData(unlist(rl, use.names=FALSE),
                                 score = unlist(lapply(resultList, score)),
                                 space = factor(names(rl)[togroup(rl)],
                                   seqlevels))
              }
              else
                gd <- do.call(c, resultList)
              if (!is.null(genome))
                genome(gd) <- genome
              gd
            } else {
              import(text = lines, format = "bed", variant = "bedGraph",
                     genome = genome, asRangedData = asRangedData,
                     which = which)
            }
        })
