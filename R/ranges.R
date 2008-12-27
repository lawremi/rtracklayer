### some utilities for working with RangesList in rtracklayer

# replace fields in 'x' with those in 'y', if given
mergeRange <- function(x, y)
{
  if (length(genome(y)))
    genome(x) <- genome(y)
  if (length(y))
    x[[1]] <- y[[1]]
  if (length(names(y)))
    names(x)[1] <- names(y)[1]
  x
}

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
