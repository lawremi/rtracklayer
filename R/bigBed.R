### =========================================================================
### BigBed support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

setClass("BigBedFile", contains = "RTLFile")
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
                           as.character(seqnames(si)), ranges(which),
                           defaultFieldIndexes, extraFieldIndexes)
            nhits <- C_ans[[1L]]
            gr <- GRanges(rep(seqnames(which), nhits), C_ans[[3L]], seqinfo=si)
            if (!is.null(C_ans[[4L]]))
              strand(gr) <- gsub(".", "*", C_ans[[4L]])
            val <- c()
            if (length(defaultFieldIndexes) && defaultFieldIndexes[1] != 0)
              val <- c(Filter(Negate(is.null), C_ans[5L:length(C_ans)]))
            val <- c(val, Filter(Negate(is.null), C_ans[[2L]]))
            elementMetadata <- DataFrame(val)
            names(elementMetadata) <- c(defaultNames ,extraNames)
            gr@elementMetadata <- elementMetadata
            gr
           })
