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
          merge(segment,
                new("genomeSegment", genome = genome, chrom = chrom,
                    start = start, end = end)))

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

# replace fields in 'x' with those in 'y', if given
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

setMethod("*", "genomeSegment",
          function(e1, e2) {
            range <- c(start(e1), end(e1))
            mid <- floor(mean(range))
            side <- diff(range)/e2/2
            start(e1) <- mid - side
            end(e1) <- mid + side
            e1
          })

setMethod("/", "genomeSegment", function(e1, e2) e1 * (1/e2))

# TODO: add intersection method

## accessors

setGeneric("genome", function(object, ...) standardGeneric("genome"))
setMethod("genome", "genomeSegment", function(object) object@genome)

setGeneric("genome<-", function(object, value) standardGeneric("genome<-"))
setReplaceMethod("genome", "genomeSegment",
                 function(object, value) {
                   object@genome <- value
                   object
                 })

setMethod("start", "genomeSegment", function(x) x@start)

setReplaceMethod("start", "genomeSegment",
           function(x, check = TRUE, value)
           {
             if (check && value > x@end)
               stop("'start' must be less than or equal to 'end'")
             x@start <- value
             x
           })

setMethod("end", "genomeSegment", function(x) x@end)

setReplaceMethod("end", "genomeSegment",
           function(x, check = TRUE, value)
           {
             if (check && x@start > value)
               stop("'start' must be less than or equal to 'end'")
             x@end <- value
             x
           })

setGeneric("chrom", function(object, ...) standardGeneric("chrom"))
setMethod("chrom", "genomeSegment", function(object) object@chrom)

setGeneric("chrom<-", function(object, value) standardGeneric("chrom<-"))
setReplaceMethod("chrom", "genomeSegment",
           function(object, value)
           {
             object@start <- value
             object
           })

## integration with Biostrings

setGeneric("genomeViews",
           function(object, segment, ...) standardGeneric("genomeViews"))
setMethod("genomeViews", c("XString", "genomeSegment"),
          function(object, segment) views(object, segment@start, segment@end))

setMethod("genomeSegment", "IRanges",
          function(object, genome = character(), chrom = character())
          genomeSegment(genome, chrom, min(start(object)), max(end(object))))
