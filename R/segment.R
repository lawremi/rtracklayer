# represents a segment of a genome
setClass("genomeSegment",
         representation(genome = "character", chrom = "character",
                        start = "numeric", end = "numeric"))

# accessors for genome segments
setGeneric("genomeSegment",
           function(object, ...) standardGeneric("genomeSegment"))
setGeneric("genomeSegment<-",
           function(object, value) standardGeneric("genomeSegment<-"))

# convenient constructors
setMethod("genomeSegment", "character",
          function(object, chrom = character(0), start = numeric(0),
                   end = numeric(0), segment = new("genomeSegment"))
          genomeSegment(genome = object, chrom = chrom, start = start,
                        end = end, segment = segment))
setMethod("genomeSegment", "missing",
          function(object, genome = character(0), chrom = character(0),
                   start = numeric(0), end = numeric(0),
                   segment = new("genomeSegment"))
          new("genomeSegment", segment, genome = genome, chrom = chrom,
              start = start, end = end))

# take the union of the segments (including gaps between segments)
c.genomeSegment <- function(...)
{
  segments <- list(...)
  stopifnot(all(lapply(segments, is, "genomeSegment")))

  first <- segments[[1]]
  if (!all(first@genome == lapply(segments, slot, "genome")) ||
      !all(first@chrom == lapply(segments, slot, "chrom")))
    warning("Genomes and/or chromosomes differ - combination not meaningful")

  genomeSegment(start = min(lapply(segments, slot, "start")),
                end = max(lapply(segments, slot, "end")), segment = first)
}

# replace empty fields in 'x' with those in 'y'
setGeneric("merge", function(x, y, ...) standardGeneric("merge"))
setMethod("merge", c("genomeSegment", "genomeSegment"),
          function(x, y)
          {
            if (length(y@genome))
              x@genome <- y@genome
            if (length(y@chrom))
              x@chrom <- y@chrom
            if (length(y@start))
              x@start <- y@start
            if (length(y@end))
              x@end <- y@end
            x
          })

# TODO: add intersection method
