setClass("ChainBlock",
         representation(ranges = "IRanges", # start in A, width
                        offset = "integer", # offset to start in B
                        score = "integer", # rle scores
                        space = "character", # rle spaces
                        reversed = "logical", # rle reversal
                        length = "integer")) # lengths for rle slots

##setMethod("length", "ChainBlock", function(x) length(x@offset))

setClass("Chain",
         prototype = prototype(elementType = "ChainBlock"),
         contains = "SimpleList")

setGeneric("import.chain",
           function(con, exclude = "_") standardGeneric("import.chain"),
           signature = "con")

setMethod("import.chain", "character", function(con, exclude) {
  .Call("readChain", con, as.character(exclude), PACKAGE="rtracklayer")
})
  
setMethod("ranges", "ChainBlock", function(x) x@ranges)
setMethod("offset", "ChainBlock", function(object) x@offset)
setMethod("score", "ChainBlock", function(x) Rle(x@score, x@length))
setMethod("space", "ChainBlock", function(x) Rle(x@space, x@length))

setGeneric("reversed", function(x, ...) standardGeneric("reversed"))
setMethod("reversed", "ChainBlock", function(x) Rle(x@reversed, x@length))

setGeneric("liftOver", function(x, chain, ...) standardGeneric("liftOver"))
setMethod("liftOver", c("GRanges", "Chain"),
          function(x, chain)
          {
            liftOverSpace <- function(ranges, chain) {
              ol <- findOverlaps(ranges, ranges(chain))
              shits <- subjectHits(ol)
              ranges <- ranges(ol, ranges, ranges(chain))
              starts <- ifelse(reversed(chain)[shits],
                               start(reflect(ranges, ranges(chain)[shits])),
                               start(ranges))
              ranges <- IRanges(starts, width=width(ranges))
              offsets <- offset(chain)[shits]
              spaces <- space(chain)[shits]
              GRanges(spaces,
                      IRanges(start(ranges) - offsets, end(ranges) - offsets))
            }
            do.call(c, mapply(liftOverSpace, x, chain))
          })
