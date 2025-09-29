### =========================================================================
### BigBed support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BigBedFile", contains = "BiocFile")
setClass("BBFile", contains = "BigBedFile")

BigBedFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("BigBedFile", resource = path)
}
BBFile <- BigBedFile

setMethod("seqinfo", "BigBedFile", function(x) {
  seqlengths <- .Call(BBDFile_seqlengths, expandPath(path(x)))
  Seqinfo(names(seqlengths), seqlengths)
})

.defaultColNames <- c("chrom", "chromStart", "chromEnd", "name", "score",
                      "strand", "thickStart", "thickEnd", "itemRgb",
                      "blockCount", "blockSizes", "blockStarts")

setClass("BigBedSelection", prototype = prototype(colnames = .defaultColNames),
         contains = "RangedSelection")

BigBedSelection <- function(ranges=IRangesList(), colnames = .defaultColNames) {
  if (is.character(ranges))
    new("BigBedSelection", GenomicSelection(ranges, colnames = colnames))
  else {
    if (is(ranges, "BigBedFile"))
      ranges <- seqinfo(ranges)
    new("BigBedSelection", ranges = as(ranges, "IntegerRangesList"),
     colnames = colnames)
  }
}

setAs("IntegerRangesList", "BigBedSelection", function(from) {
  new("BigBedSelection", as(from, "RangedSelection"), colnames = .defaultColNames)
})

setAs("GenomicRanges", "BigBedSelection", function(from) {
  as(as(from, "IntegerRangesList"), "BigBedSelection")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.bb", function(con, ...) standardGeneric("import.bb"))

setMethod("import.bb", "ANY", function(con, ...) {
  import(con, "BigBed", ...)
})

setMethod("import", "BigBedFile",
          function(con, format, text, selection = BigBedSelection(which, ...),
                   which = con, ...) {

            if (!missing(format))
              checkArgFormat(con, format)
            si <- seqinfo(con)

            # If which is NULL, create a GRanges spanning all sequences full length
            if (is.null(which)) {
              full_ranges <- IRanges(start = rep(1, length(seqlengths(si))),
                                    end = seqlengths(si))
              which <- GRanges(seqnames = seqlevels(si), ranges = full_ranges)
              selection <- BigBedSelection(which)
            } else {
              selection <- as(selection, "BigBedSelection")
            }

            col_names <- colnames(selection)
            # Get all available field names from bigBed file
            fields_name <- .Call(BBDFile_fieldnames, expandPath(path(con)))
            fields_name <- c(fields_name[[1L]], fields_name[[2L]])
            # always prefer file field names if the user did not provide any
            if (identical(col_names, .defaultColNames))
                col_names <- fields_name
            required_fields <- c("chrom", "chromStart", "chromEnd")
            missing_fields <- setdiff(required_fields, col_names)
            col_names <- c(col_names, missing_fields)
            # set the selected colnames back
            colnames(selection) <- col_names
            unavailable_fields <- setdiff(col_names, fields_name)
            unavailable_fields <- setdiff(unavailable_fields, .defaultColNames)
            if (length(unavailable_fields) > 0L)
                warning(paste("Unavailable field(s):",
                        paste(unavailable_fields, collapse = ", ")))

            # Extract and sanitize ranges against seqlevels in file
            ranges <- ranges(selection)
            badSpaces <- setdiff(names(ranges)[lengths(ranges) > 0L], seqlevels(si))
            if (length(badSpaces) > 0L)
              warning("'which' contains seqnames not known to BigBed file: ",
                      paste(badSpaces, collapse = ", "))
            ranges <- ranges[names(ranges) %in% seqlevels(si)]

            # Flatten and split by seqlevels for querying
            flatranges <- unlist(ranges, use.names = FALSE)
            if (is.null(flatranges))
              flatranges <- IRanges()
            which_rl <- split(flatranges, factor(space(ranges), seqlevels(si)))
            which <- GRanges(which_rl)

            C_ans <- .Call(BBDFile_query, expandPath(path(con)),
                           as.character(seqnames(which)), ranges(which),
                           colnames(selection))


            nhits <- C_ans[["n_qhits"]]
            # Reconstruct GRanges with genomic coordinates and seqinfo
            result_seqnames <- C_ans[["chrom"]]
            chromStart <- C_ans[["chromStart"]]
            chromEnd <- C_ans[["chromEnd"]]
            # to 1 based
            chromStart <- chromStart + 1L
            result_ranges <- IRanges(start = chromStart, chromEnd)
            gr <- GRanges(result_seqnames, result_ranges, seqinfo = si)

            if ("strand" %in% names(C_ans) && !is.null(C_ans[["strand"]])) {
                cleaned_strand <- gsub("\\.", "*", C_ans[["strand"]])
                strand(gr) <- factor(cleaned_strand, levels = c("+", "-", "*"))
            }

            if ("thickStart" %in% names(C_ans) && !is.null(C_ans[["thickStart"]])) {
                # convert 0-based to 1-based start
                thickStart <- C_ans[["thickStart"]] + 1L
                thickEnd <- C_ans[["thickEnd"]]

                if (!is.null(thickEnd) && length(thickStart) == length(thickEnd)) {
                    thick <- IRanges(start = thickStart, end = thickEnd)
                    mcols(gr)$thick <- thick
                } else {
                    warning("thickEnd not found or length mismatch with thickStart")
                }
            }

            if ("blockCount" %in% names(C_ans) && !is.null(C_ans[["blockCount"]])) {
                blockCount <- C_ans[["blockCount"]]
                blockSizes <- C_ans[["blockSizes"]]
                blockStarts <- C_ans[["blockStarts"]]
                # Convert blockStarts from 0-based to 1-based coordinates
                ir_list <- IRanges(start = blockStarts + 1L, width = blockSizes)
                blocks <- relist(ir_list, PartitioningByWidth(blockCount))
                names(blocks) <- NULL
                mcols(gr)$blocks <- blocks
            }

            if ("itemRgb" %in% names(C_ans) && !is.null(C_ans[["itemRgb"]])) {
                color <- C_ans[["itemRgb"]]
                spec <- color != "0"
                cols <- unlist(strsplit(color[spec], ",", fixed=TRUE),
                               use.names=FALSE)
                cols <- matrix(as.integer(cols), 3)
                color <- rep(NA, length(gr))
                color[spec] <- rgb(cols[1,], cols[2,], cols[3,],
                                   maxColorValue = 255L)
                mcols(gr)$itemRgb <- color
            }

            processed_fields <- c("n_qhits", "chrom", "chromStart", "chromEnd",
                                  "strand", "thickStart", "thickEnd",
                                  "blockCount", "blockSizes", "blockStarts",
                                  "itemRgb")
            remaining_fields <- setdiff(names(C_ans), processed_fields)
            for (field in remaining_fields)
                mcols(gr)[[field]] <- C_ans[[field]]

           gr
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.bb", function(object, con, ...) standardGeneric("export.bb"))

setMethod("export.bb", "ANY",
          function(object, con, ...)
          {
            export(object, con, "BigBed", ...)
          })

setMethod("export", c("ANY", "BigBedFile"),
          function(object, con, format, ...)
          {
            object <- as(object, "GRanges")
            callGeneric()
          })

setMethod("export", c("GenomicRanges", "BigBedFile"),
          function(object, con, format, compress = TRUE, extraIndexes = "")
          {
            if (!missing(format))
              checkArgFormat(con, format)
            con <- path.expand(path(con))
            object <- sortBySeqnameAndStart(object)
            seqlengths <- seqlengths(object)
            stopIfNotValidForExport(object)
            if (!is.character(extraIndexes))
              stop("The extraIndexes must be character")
            if (any(is.na(seqlengths)))
              stop("Unable to determine seqlengths; either specify ",
                   "'seqlengths' or specify a genome on 'object' that ",
                   "is known to BSgenome or UCSC")
            if (!isTRUEorFALSE(compress))
              stop("'compress' must be TRUE or FALSE")
            seqlengths <- seqlengths(object)
            bedString <- bedString(object)
            autoSqlString <- autoSqlString(object)
            extraIndexes <- gsub("[\n\t ]", "", extraIndexes, perl = TRUE)
            invisible(BigBedFile(.Call(BBDFile_write, seqlengths, bedString, autoSqlString,
                                       extraIndexes, compress, con)))
          })

stopIfNotValidForExport <- function(x) {
  elementMetadata <- elementMetadata(x)
  name <- elementMetadata$name
  score <- elementMetadata$score
  itemRgb <- elementMetadata$itemRgb
  thick <- elementMetadata$thick
  blocks <- elementMetadata$blocks
  if (!is.null(name) && (!is.character(name) || any(is.na(name))))
    stop("The name must be character, without any NA's")
  if (isValidScore(score))
    stop("The score must be numeric, without any NA's")
  if (!is.null(itemRgb) && (!is.character(itemRgb) || any(is.na(itemRgb))))
    stop("The itemRgb must be character, without any NA's")
  if (!is.null(thick) && !is(thick, "IRanges"))
    stop("The thick must be IRanges")
  if (!is.null(blocks) && !is(blocks, "IRangesList"))
    stop("The blocks must be IRangesList")
}

bedString <- function(x) {
  elementMetadata <- elementMetadata(x)
  name <- elementMetadata$name
  elementMetadata$name <- NULL
  score <- elementMetadata$score
  elementMetadata$score <- NULL
  strand <- as.character(strand(x))
  strand <- gsub("*", ".", strand, fixed = TRUE)
  thick <- elementMetadata$thick
  thickStart <- NULL
  thickEnd <- NULL
  if (!is.null(thick)) {
    thickStart <- start(ranges(thick)) - 1L
    thickEnd <- end(ranges(thick))
    elementMetadata$thick <- NULL
  }
  itemRgb <- as.data.frame(t(col2rgb(elementMetadata$itemRgb)))
  itemRgb <- do.call(paste, c(itemRgb, sep=","))
  elementMetadata$itemRgb <- NULL
  blocks <- elementMetadata$blocks
  blockCount  <- NULL
  blockSizes  <- NULL
  blockStarts <- NULL
  if (!is.null(blocks)) {
    length <- length(blocks)
    blockCount <- lengths(blocks)
    blockSizes <- lapply(width(blocks), function(x) paste(x, collapse=","))
    blockStarts <- lapply(start(blocks) - 1L, function(x) paste(x, collapse=","))
    elementMetadata$blocks <- NULL
  }
  extraColumnsString <- do.call(paste, as.list(elementMetadata))
  paste(as.character(seqnames(x)), start(ranges(x)) - 1L, end(ranges(x)), name, score,
                     strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes,
                     blockStarts, extraColumnsString, collapse = "\n")
}

autoSqlString <- function(x) {
  asString <- c('table bed "Browser Extensible Data" (\n',
                'string chrom; "Reference sequence chromosome or scaffold"\n',
                'uint chromStart; "Start position in chromosome"\n',
                'uint chromEnd; "End position in chromosome"\n')

  names <- c("name", "itemRgb", "score", "thick", "blocks", "double", "integer", "character", "raw")
  values <- c('string name; "Name of item."\n',
              'uint reserved; "Used as itemRgb as of 2004-11-22"\n',
              'uint score; "Score (0-1000)"\nchar[1] strand; "+ or - for strand"\n',
              paste0('uint thickStart; "Start of where display should be thick (start codon)"\n',
                     'uint thickEnd; "End of where display should be thick (stop codon)"\n'),
              paste0('int blockCount; "Number of blocks"\n',
                     'int[blockCount] blockSizes; "Comma separated list of block sizes"\n',
                     'int[blockCount] chromStarts; "Start positions relative to chromStart"\n'),
              "double ", "int ", "string ", "uint ")
  mapping <- setNames(values, names)
  metadata <- elementMetadata(x)
  names <- names(metadata)
  defaultFields <- colnames(BigBedSelection())
  fieldsString <- lapply(names, function(y) {
    if (y %in% defaultFields)
      mapping[y]
    else {
      typeString <- mapping[storage.mode(metadata[[y]])]
      paste(typeString, y, '; ""\n')
    }
  })
  asString <- c(asString, fieldsString, ')')
  paste(asString, collapse = "")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

cleanupBigBedCache <- cleanupBigWigCache
