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

.defaultColNames <- c("name", "score", "thick", "itemRgb", "blocks")

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
                   which = con, ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            si <- seqinfo(con)
            selection <- as(selection, "BigBedSelection")
            ranges <- ranges(selection)
            badSpaces <- setdiff(names(ranges)[lengths(ranges) > 0L], seqlevels(si))
            if (length(badSpaces) > 0L)
              warning("'which' contains seqnames not known to BigBed file: ",
                      paste(badSpaces, collapse = ", "))
            ranges <- ranges[names(ranges) %in% seqlevels(si)]
            flatranges <- unlist(ranges, use.names = FALSE)
            if (is.null(flatranges))
              flatranges <- IRanges()
            which_rl <- split(flatranges, factor(space(ranges), seqlevels(si)))
            which <- GRanges(which_rl)
            allFields <- .Call(BBDFile_fieldnames, expandPath(path(con)))
            defaultFields <- allFields[[1L]]
            ValidextraFields <- allFields[[2L]]
            selectedFields <- colnames(selection)
            extraFields <- setdiff(selectedFields, defaultFields)
            if (identical(colnames(BigBedSelection()), selectedFields)) {
              selectedFields <- defaultFields[defaultFields != ""]
              extraFields <- ValidextraFields[ValidextraFields != ""]
            }
            if (!identical(selectedFields, defaultFields)) {
              defaultFields <- defaultFields[defaultFields != ""]
              ValidextraFields <- ValidextraFields[ValidextraFields != ""]
              defaultFieldIndexes <- which(defaultFields %in% selectedFields)
              extraFieldIndexes <- which(ValidextraFields %in% extraFields)
              invalidFields <- setdiff(extraFields, ValidextraFields)
              if (length(defaultFieldIndexes) == 0L)
                defaultFieldIndexes <- c(0L)
              if (length(extraFieldIndexes) == 0L)
                extraFieldIndexes <- c(0L)
              if (length(invalidFields))
                warning("Invalid ", invalidFields, " field(s)")
            }else {
              defaultFieldIndexes <- c()
              extraFieldIndexes <- c()
            }
            defaultNames <- defaultFields[defaultFields %in% selectedFields]
            extraNames <- ValidextraFields[ValidextraFields %in% extraFields]
            C_ans <- .Call(BBDFile_query, expandPath(path(con)),
                           as.character(seqnames(which)), ranges(which),
                           defaultFieldIndexes, extraFieldIndexes)
            nhits <- C_ans[[1L]]
            gr <- GRanges(rep(seqnames(which), nhits), C_ans[[3L]], seqinfo=si)
            if (!is.null(C_ans[[4L]]))
              strand(gr) <- gsub(".", "*", C_ans[[4L]], fixed = TRUE)
            blocksPosition <- which(defaultNames %in% c("blocks"))
            if (length(blocksPosition)) {
              blocksPosition <- 4 + blocksPosition
              C_ans[[blocksPosition]] <- IRangesList(C_ans[[blocksPosition]])
            }
            val <- c()
            if (length(defaultFieldIndexes) && defaultFieldIndexes[1] != 0)
              val <- c(Filter(Negate(is.null), C_ans[5L:length(C_ans)]))
            val <- c(val, Filter(Negate(is.null), C_ans[[2L]]))
            elementMetadata <- DataFrame(val)
            names(elementMetadata) <- c(defaultNames ,extraNames)
            gr@elementMetadata <- elementMetadata
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
    thickStart <- start(ranges(thick))
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
    blockStarts <- lapply(start(blocks), function(x) paste(x, collapse=","))
    elementMetadata$blocks <- NULL
  }
  extraColumnsString <- do.call(paste, as.list(elementMetadata))
  paste(as.character(seqnames(x)), start(ranges(x)), end(ranges(x)), name, score,
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
