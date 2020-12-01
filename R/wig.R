### =========================================================================
### WIG (Wiggle) support (fixedStep and variableStep, legacy bedGraph)
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("WIGFile", contains = "BiocFile")

WIGFile <- function(resource) {
  new("WIGFile", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.wig",
           function(object, con, ...) standardGeneric("export.wig"))

setMethod("export.wig", "ANY",
          function(object, con, ...)
          {
            export(object, con, "wig", ...)
          })

setMethod("export", c("ANY", "WIGFile"),
          function(object, con, ...)
          {
            cl <- class(object)
            track <- try(as(object, "GRanges"), silent = TRUE)
            if (class(track) == "try-error") {
              track <- try(as(object, "SimpleGRangesList"), silent = TRUE)
              if (is(track, "try-error"))
                stop("cannot export object of class '", cl, "': ", track)
            }
            export(track, con, ...)
          })

.wigWriter <- function(chromData, con, dataFormat, append) {
  m <- manager()
  con <- connection(m, con, if (append) "a" else "w")
  on.exit(release(m, con))
  cat(dataFormat, file = con)
  cat(" chrom=", as.character(seqnames(chromData)[1]), file = con, sep = "")
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

isValidScore <- function(score) {
    !(is.numeric(score) || is(score, "Rle") && is.numeric(runValue(score))) ||
        any(is.na(score))
}

setMethod("export", c("GenomicRanges", "WIGFile"),
          function(object, con, format,
                   dataFormat = c("auto", "variableStep", "fixedStep"),
                   writer = .wigWriter, append = FALSE, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (isValidScore(score(object)))
              stop("The score must be numeric, without any NA's")
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            dataFormat <- match.arg(dataFormat)            
            doBlock <- function(chromData) {
              if (length(chromData) == 0L)
                  return()
              if (is.unsorted(start(chromData)))
                chromData <- chromData[order(start(chromData)),]
              starts <- start(chromData)
              ends <- end(chromData)
              if (!all(tail(starts, -1) - head(ends, -1) > 0))
                stop("Features cannot overlap. ",
                     "Note that WIG does not distinguish between strands - ",
                     "try exporting two tracks, one for each strand.")
              if (length(starts) == 1)
                steps <- 0
              else steps <- diff(starts)
              fixedSpan <- is(object, "GPos") || {
                  spans <- ends - starts + 1L
                  all(spans[1L] == spans)
              }
              if (!fixedSpan)
                stop("The span must be uniform for Wiggle export. ",
                     "Consider exporting to bedGraph/bigWig, ",
                     "or coerce data to a GPos object first.")
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
            lapply(split(object, seqnames(object)), doBlock)
            invisible(con)
          })

setMethod("export", c("UCSCData", "WIGFile"),
          function(object, con, format, trackLine = TRUE, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (trackLine) {
              export.ucsc(object, con, ...)
            } else {
              callNextMethod()
            }
            invisible(con)
          })

setMethod("export", c("SimpleGRangesList", "WIGFile"),
          .export_SimpleGRangesList_BiocFile)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.wig", function(con, ...) standardGeneric("import.wig"))

setMethod("import.wig", "ANY",
          function(con, ...)
          {
            import(con, "wig", ...)
          })

setMethod("import", "WIGFile",
          function(con, format, text, genome = NA,
                   trackLine = TRUE, which = NULL, seqinfo = NULL, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            file <- con
            m <- manager()
            con <- connection(m, con, "r")
            on.exit(release(m, con))
            ## check for a track line
            line <- scanTrackLine(con)
            if (!is.null(line) && trackLine) {
              pushBack(line, con)
              return(import.ucsc(initialize(file, resource = con), drop = TRUE,
                                 trackLine = FALSE, genome = genome,
                                 which = which, seqinfo = seqinfo,
                                 ...))
            }
            lines <- readLines(con, warn = FALSE)
            formatInds <- grep("chrom=", lines, fixed=TRUE)
            formatLines <- lines[formatInds]
            starts <- formatInds + 1L
            ends <- c(tail(formatInds, -1) - 1L, length(lines))
            format <- gsub("^([^ ]*) .*", "\\1", formatLines)
            parsedFormat <- lapply(formatLines, ucscParsePairs)
            chrom <- sapply(parsedFormat, `[[`, "chrom")
            if (is.null(seqinfo) && !is.null(genome) && !is.na(genome))
              seqinfo <- seqinfoForGenome(genome)
            if (!is.null(seqinfo))
              seqlevels <- seqlevels(seqinfo)
            else seqlevels <- unique(chrom)
            if (length(formatLines)) {
              parseData <- function(i) {
                ## parse format line
                formatVals <- parsedFormat[[i]]
                # parse the data values
                block_lines <- window(lines, starts[i], ends[i])
                if (!length(block_lines)) {
                  gr <- GRanges(score = numeric())
                  seqlevels(gr) <- seqlevels
                  return(gr)
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
                                  which = which)
                seqlevels(gd) <- seqlevels
                gd
              }
              parseInds <- seq_along(formatInds)
              if (!is.null(which)) {
                message("For efficiency, consider converting this WIG file",
                        " to a BigWig file;\nsee ?wigToBigWig")
                parseInds <- parseInds[chrom %in% seqlevels(which)]
              }
              resultList <- lapply(parseInds, parseData)
              gd <- do.call(c, resultList)
              if (!is.null(seqinfo))
                seqinfo(gd) <- seqinfo
              else if (!is.null(genome))
                genome(gd) <- genome
              gd
            } else {
              import(text = lines, format = "bedGraph",
                     genome = genome, which = which, seqinfo = seqinfo)
            }
        })
