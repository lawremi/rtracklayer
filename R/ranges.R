### some utilities for working with RangesList in rtracklayer

# replace fields in 'x' with those in 'y', if given
setGeneric("merge", function(x, y, ...) standardGeneric("merge"))
setMethod("merge", c("RangesList", "RangesList"),
          function(x, y)
          {
            if (length(genome(y)))
              genome(x) <- genome(y)
            if (length(names(y)))
              names(x)[1] <- names(y)[1]
            x[[1]] <- y[[1]]
            x
          })

setMethod("*", "RangesList",
          function(e1, e2) {
            r <- e1[[1]]
            if (e2 < 0)
              e2 <- abs(1/e2)
            range <- c(start(r), end(r))
            mid <- floor(mean(range))
            side <- (diff(range)+1)/e2/2
            start(r) <- mid - side
            end(r) <- mid + side
            e1[[1]] <- r
            e1
          })

### Is this used? With Ranges, can call range(c(...))
## take the union of the segments (including gaps between segments)
## c.genomeSegment <- function(...)
## {
##   segments <- list(...)
##   stopifnot(all(lapply(segments, is, "genomeSegment")))

##   first <- segments[[1]]
##   if (!all(first@genome == lapply(segments, slot, "genome")) ||
##       !all(first@chrom == lapply(segments, slot, "chrom")))
##     warning("Genomes and/or chromosomes differ - combination not meaningful")

##   genomeSegment(start = min(lapply(segments, slot, "start")),
##                 end = max(lapply(segments, slot, "end")), segment = first)
## }
