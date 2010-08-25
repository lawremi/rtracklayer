setClass("AlignmentSpace",
         representation(ranges = "IRanges", # start in A, width
                        offset = "integer", # offset to start in B
                        score = "integer", # rle scores
                        space = "character", # rle spaces
                        rev = "logical", # rle reversal
                        length = "integer")) # lengths for rle slots

##setMethod("length", "AlignmentSpace", function(x) length(x@offset))

setClass("Alignment",
         prototype = prototype(elementType = "AlignmentSpace"),
         contains = "SimpleList")

read.chain <- function(path, exclude = "_") {
  .Call("readChain", path, exclude, PACKAGE="IRanges")
}

setMethod("score", "AlignmentSpace", function(x) x@score)

setGeneric("map", function(x, alignment, ...) standardGeneric("map"))
setMethod("map", c("RangesList", "Alignment"),
          function(x, alignment)
          {
            r <- IRanges()
            s <- character()
            for (space in names(x)) {
              ranges <- x[[space]]
              align <- alignment[[space]]
              ol <- findOverlaps(ranges, ranges(align))
              hits <- as.matrix(ol)
              ranges <- ranges(ol, ranges, ranges(align))
              starts <- ifelse(reversed(align)[hits[,2L]],
                               start(reflect(ranges, ranges(align)[hits[,2L]])),
                               start(ranges))
              ranges <- IRanges(starts, width=width(ranges))
              offsets <- offset(align)[hits[,2L]]
              spaces <- space(align)[hits[,2L]]
              r <- c(r, IRanges(start(ranges) - offsets, end(ranges) - offsets))
              s <- c(s, spaces)
            } ### FIXME: need some more efficient way of bundling result
            split(r, s)
          })
