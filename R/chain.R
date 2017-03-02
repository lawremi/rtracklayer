### =========================================================================
### Chain file parsing and lift over
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes
###

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

setClass("ChainFile", contains = "RTLFile")

ChainFile <- function(path) {
  if (!isSingleString(path))
    stop("'filename' must be a single string, specifying a path")
  new("ChainFile", resource = path)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.chain",
           function(con, ...) standardGeneric("import.chain"))

setMethod("import.chain", "ANY", function(con, ...) {
  import(con, "chain", ...)
})

setMethod("import", "ChainFile", function(con, format, text, exclude = "_") {
  if (!missing(format))
    checkArgFormat(con, format)
  if (!is(connection(con), "file")) {
    stop("chain import currently only handles local, uncompressed files")
  }
  .Call("readChain", path.expand(path(con)), as.character(exclude),
        PACKAGE="rtracklayer")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("ranges", "ChainBlock", function(x) x@ranges)
setMethod("offset", "ChainBlock", function(object) object@offset)
setMethod("score", "ChainBlock", function(x) Rle(x@score, x@length))
setMethod("space", "ChainBlock", function(x) Rle(x@space, x@length))

setGeneric("reversed", function(x, ...) standardGeneric("reversed"))
setMethod("reversed", "ChainBlock", function(x) Rle(x@reversed, x@length))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Liftover
###

flipStrandSimple <- function(strand, flip) {
  strand <- as.vector(strand)
  flipped <- ifelse(flip, ifelse(strand == "+", "-",
                                 ifelse(strand == "-", "+", strand)),
                    strand)
  strand(flipped)
}

flipStrandTricky <- function(strand, flip) {
  strandCodes <- c("+" = 1L, "-" = -1L, "*" = 0L)
  strandInt <- strandCodes[as.vector(strand)]
  flipped <- ifelse(flip, strandInt * -1L, strandInt) + 2L
  strandRevCodes <- factor(c("-", "*", "+"), levels(strand()))
  strandRevCodes[as.vector(flipped)]
}

setGeneric("liftOver", function(x, chain, ...) standardGeneric("liftOver"))
setMethod("liftOver", c("GenomicRanges", "Chain"),
          function(x, chain)
          {
            liftOverSpace <- function(gr, chain, subind) {
              r <- ranges(gr)
              ol <- findOverlaps(r, ranges(chain))
              shits <- subjectHits(ol)
              r <- ranges(ol, r, ranges(chain))
              rev <- as.vector(reversed(chain)[shits])
              starts <- ifelse(rev,
                               start(reflect(r, ranges(chain)[shits])),
                               start(r))
              strand <- flipStrandTricky(strand(gr)[queryHits(ol)], rev)
              r <- IRanges(starts, width=width(r))
              offsets <- offset(chain)[shits]
              spaces <- space(chain)[shits]
              ind[[as.character(seqnames(gr)[1])]] <<- subind[queryHits(ol)]
              GRanges(spaces,
                      IRanges(start(r) - offsets, end(r) - offsets),
                      strand = strand,
                      values(gr)[queryHits(ol),,drop=FALSE])
            }
            rl <- split(x, seqnames(x), drop = TRUE)
            unchainedNames <- setdiff(names(rl), names(chain))
            if (length(unchainedNames))
              message("Discarding unchained sequences: ",
                      paste(unchainedNames, collapse = ", "))
            sharedNames <- intersect(names(rl), names(chain))
            ind <- split(seq_len(length(x)),
                         as.vector(seqnames(x)))[sharedNames]
            liftedList <- mapply(liftOverSpace, rl[sharedNames],
                                 chain[sharedNames], ind, SIMPLIFY=FALSE)
            lifted <- unlist(GRangesList(liftedList), use.names=FALSE)
            f <- structure(unlist(ind, use.names=FALSE),
                           levels=seq_len(length(x)), class="factor")
            setNames(split(lifted, f), names(x))
          })

setMethod("liftOver", c("ANY", "ANY"),
          function(x, chain) {
    chain <- as(chain, "Chain")
    callGeneric()
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###


## mismatches take into account dels but not start, need opposite
mapMismatches <- function(mismatch, range, dels) {
    names(mismatch) <- seq_along(mismatch)
    mmflat <- unlist(mismatch)
    mmgr <- GRanges(names(mmflat), IRanges(mmflat, width=1L))
    seqlevels(mmgr) <- names(mismatch)
    start(mapToUnaligned(mmgr, range, dels))
}

mapToUnaligned <- function(x, range, dels) {
    tx <- gaps(dels, 1L, end(range) + sum(width(dels)))
    txgr <- as(setNames(tx, seqlevels(x)), "GRanges")
    txgrl <- split(txgr, seqnames(txgr))
    mmlocal <- mapToTranscripts(x, txgrl)
    ##seqlevels(mmlocal) <- seqlevels(x)
    shift(split(mmlocal, seqnames(mmlocal)), start(range) - 1L)
}

deletionsFromGaps <- function(x) {
    gapraw <- as.raw(DNAString("-"))
    as(RleList(as(x, "RawList") == gapraw), "IRangesList")
}

mapContiguousRanges <- function(x, aln) {
    ranges(unlist(mapToUnaligned(x, ranges(aln), indel(aln)), use.names=FALSE))
}

setAs("AlignedXStringSetList", "Chain", function(from) {
    pairs <- zipdown(from[lengths(from) == 2L])
    dels <- c(indel(first(pairs)), indel(second(pairs)))
    names(dels) <- rep(seq_along(pairs), 2)
    delgr <- as(dels, "GRanges")
    seqlengths(delgr) <-
        width(first(pairs)) + sum(width(indel(first(pairs))))          
    dr <- disjoin(subset(gaps(delgr), strand=="*"))
    r1 <- mapContiguousRanges(dr, first(pairs))
    r2 <- mapContiguousRanges(dr, second(pairs))
    block <- new("ChainBlock",
                 ranges = r1,
                 offset = start(r1) - start(r2),
                 score = rep(NA_integer_, length(pairs)),
                 space = rep("second", length(pairs)),
                 reversed =
                     inverted(first(pairs)) != inverted(second(pairs)),
                 length = as.integer(table(seqnames(dr))))
    new("Chain", SimpleList(first=block))
})

