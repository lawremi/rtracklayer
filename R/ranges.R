### =========================================================================
### Genome-oriented methods for RangedData classes
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("genome", function(x, ...) standardGeneric("genome"))
setMethod("genome", "RangedData", function(x) universe(x))
setMethod("genome", "GRanges",
          function(x) {
            if (is.null(metadata(x)) || is.character(metadata(x))) 
              metadata(x)
            else
              metadata(x)$genome
          })

setGeneric("genome<-", function(x, value) standardGeneric("genome<-"))
setReplaceMethod("genome", "RangedData",
                 function(x, value) {
                   genome(x@ranges) <- value
                   x
                 })
setReplaceMethod("genome", "GRanges",
                 function(x, value) {
                   if (!is.null(value) && !IRanges:::isSingleString(value)) 
                     stop("'value' must be a single string or NULL")
                   metadata(x)$genome <- value
                   x
                 })

setGeneric("chrom", function(x, ...) standardGeneric("chrom"))
setMethod("chrom", "RangedData", function(x) chrom(ranges(x)))
setMethod("chrom", "GRanges", function(x) seqnames(x))

setGeneric("chrom<-", function(x, ..., value) standardGeneric("chrom<-"))
setReplaceMethod("chrom", "RangedData", function(x, value) {
  chrom(ranges(x)) <- value
  x
})
setReplaceMethod("chrom", "GRanges", function(x, value) {
  seqnames(x) <- value
  x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicData <- function(ranges, ..., strand = NULL, chrom = NULL, genome = NULL,
                        asRangedData = TRUE)
{
  if (!IRanges:::isTRUEorFALSE(asRangedData))
    stop("'asRangedData' must be TRUE or FALSE")
  if (!is(ranges, "Ranges")) {
    if (is(ranges, "data.frame") || is(ranges, "DataTable")) {
      colnames(ranges)[match("chrom", colnames(ranges))] <- "space"
    }
    gd <- RangedData(ranges) # direct coercion
    if (!asRangedData)
      gd <- as(gd, "GRanges")
  } else {
    if (length(chrom) > length(ranges))
      stop("length of 'chrom' greater than length of 'ranges'")
    if (length(chrom) > 0 && (length(ranges) %% length(chrom) != 0))
      stop("length of 'ranges' not a multiple of 'chrom' length")
    if (!is.null(genome) && (length(genome) != 1 || !is.character(genome)))
      stop("'genome' must be a single string")
    if (asRangedData) {
      if (!is.null(strand)) {
        if (!all(strand[!is.na(strand)] %in% levels(strand())))
          stop("strand values should be 'NA', '-', '+' or '*'")
        strand <- factor(strand, levels(strand()))
        gd <- RangedData(ranges, ..., strand = strand, space = chrom,
                         universe = genome)
      } else {
        gd <- RangedData(ranges, ..., space = chrom, universe = genome)
      }
    } else {
      if (is.null(chrom))
        chrom <- Rle(factor("1"), length(ranges))
      dots <- list(...)
      if (length(dots) == 1) {
        dots <- dots[[1L]]
        if ((is(dots, "data.frame") || is(dots, "DataTable")) &&
            !is.null(dots[["strand"]])) {
          strand <- dots[["strand"]]
          dots[["strand"]] <- NULL
          return(GenomicData(ranges = ranges, dots, strand = strand,
                             chrom = chrom, genome = genome,
                             asRangedData = asRangedData))
        }
      }
      if (is.null(strand))
        strand <- Rle("*", length(ranges))
      gd <- GRanges(seqnames = chrom, ranges = ranges, strand = strand, ...)
      if (!is.null(genome))
        metadata(gd) <- list(universe = genome)
    }
  }
  gd
}

### =========================================================================
### Genome-oriented methods for GRanges/RangesList classes
### -------------------------------------------------------------------------

setMethod("genome", "RangesList", function(x) universe(x))

setReplaceMethod("genome", "RangesList",
                 function(x, value) {
                   universe(x) <- value
                   x
                 })

setMethod("chrom", "RangesList", function(x) {
  chrom <- names(x)
  if (!is.null(chrom))
    chrom <- rep(factor(chrom, chrom), unlist(lapply(x, length)))
  chrom
})

setReplaceMethod("chrom", "RangesList", function(x, value) {
  names(x) <- value
  x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### DEPRECATED
GenomicRanges <- function(start = integer(), end = integer(), chrom = NULL,
                          genome = NULL)
{
  .Deprecated("GRangesForBSGenome or GRangesForUCSCGenome")
  ir <- IRanges(start, end)
  if (!is.null(chrom)) {
    if (!is.factor(chrom))
      chrom <- factor(chrom, unique(chrom))
    if (!length(ir))
      chrom <- chrom[FALSE]
    if (length(chrom) != length(ir))
      stop("length of 'chrom' must equal that of 'start' and 'end'",
           " unless 'start' and 'end' are zero length")
    rl <- split(ir, chrom)
  } else rl <- RangesList(ir)
  universe(rl) <- genome
  rl
}

GRangesForBSGenome <- function(genome, chrom = NULL, start = 1L, end = NULL,
                               ...)
{
  if (missing(genome) || !IRanges:::isSingleString(genome))
    stop("'genome' must be a single string identifying a genome")
  bsgenome <- .genomeForID(genome)
  if (is.null(bsgenome))
    stop("genome '", genome,
         "' does not correspond to an installed BSgenome package")
  GRangesForGenome(genome, seqlengths(bsgenome), chrom = chrom, start = start,
                   end = end, ...)
}

## Internal helper
GRangesForGenome <- function(genome, seqlens, chrom = NULL, start = 1L,
                             end = NULL, ...)
{
  if (!is.integer(seqlens) || is.null(names(seqlens)))
    stop("'seqlens' must be a named integer vector of chromosome lengths")
  if (any(start < 1))
    stop("all values in 'start' should be positive")
  if (is.null(chrom))
    chrom <- names(seqlens)
  else {
    badChrom <- setdiff(chrom, names(seqlens))
    if (length(badChrom))
      stop("Chromosome(s) ", paste(badChrom, collapse = ", "),
           "are invalid for: ", genome)
  }
  if (is.null(end))
    end <- seqlens[chrom]
  gr <- GRanges(chrom, IRanges(start, end), ..., seqlengths = seqlens)
  genome(gr) <- genome
  gr
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

## normalize 'range', using 'session' for default genome
normGenomeRange <- function(range, session) {
  ## the user can specify a portion of the genome in several ways:
  ## - String identifying a genome
  ## - RangesList, possibly constructed with deprecated GenomicRanges()
  ## - GRanges, the preferred way, possibly from GRangesForUCSCGenome()
  ## - We do not allow Ranges, since it does not make sense to have one range
  ##   over many chromosomes
  if (is.character(range)) {
    genome(session) <- range
    return(GRangesForGenome(range, seqlengths(session)))
  }
  genome <- genome(session)  
  if (is.null(genome(range)))
    genome(range) <- genome
  else if (genome(range) != genome) {
    genome(session) <- genome(range)
    on.exit(genome(session) <- genome)
  }
  if (is(range, "RangesList")) {
    warning("Specifying coordinates with 'RangesList' is deprecated: ",
            "Use 'GRanges' instead")
    chrom <- names(range)
    start <- 1L
    end <- NULL
    if (!is.null(chrom)) {
      flatRange <- unlist(range)
      if (length(flatRange)) {
        start <- start(flatRange)
        end <- end(flatRange)
      }
    }
    GRangesForGenome(genome(range), seqlengths(session), chrom, start, end)
  } else if (is(range, "GRanges")) {
    badChroms <- setdiff(seqnames(range), seqnames(session))
    if (length(badChroms))
      stop("Invalid chromosomes for ", genome(range), ": ",
           paste(badChroms, collapse = ", "))
    range
  }
  else stop("'range' should be either a genome string, RangesList or GRanges")
}

spansGenome <- function(x) {
  strand(x) <- "*"
  x <- reduce(x)
  w <- structure(width(x), names = as.character(seqnames(x)))
  identical(w[names(seqlengths(x))], seqlengths(x))
}

### =========================================================================
### Genome-oriented conveniences for RangedSelection classes
### -------------------------------------------------------------------------

.genomeForID <- function(genome) {
  pkgs <- grep("^BSgenome\\.", rownames(installed.packages()), value = TRUE)
  pkg <- grep(paste(genome, "$", sep = ""), pkgs, value = TRUE)
  if (length(pkg) == 1) {
    org <- strsplit(pkg, ".", fixed=TRUE)[[1]][2]
    get(org, getNamespace(pkg))
  } else NULL
}

## One could imagine the BSgenome object having a coerce method to
## Ranges and RangesList, and this function could use that. But would
## the coercion consider masks?

GenomicSelection <- function(genome, chrom = NULL, colnames = character(0))
{
  if (missing(genome) || !IRanges:::isSingleString(genome))
    stop("'genome' must be a single string identifying a genome")
  bsgenome <- .genomeForID(genome)
  if (is.null(bsgenome))
    stop("genome '", genome,
         "' does not correspond to an installed BSgenome package")
  lens <- seqlengths(bsgenome)
  if (is.null(chrom))
    chrom <- names(lens)
  else {
    if (!is.character(chrom))
      stop("'chrom' must be NULL or a character vector")
    invalidChroms <- setdiff(chrom, names(lens))
    if (length(invalidChroms))
      stop("'chrom' contains invalid chromosomes: ",
           paste(invalidChroms, collapse = ", "))
    lens <- lens[chrom]
  }
  RangedSelection(split(IRanges(1, lens), factor(chrom, chrom)), colnames)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Questionable utilities
###

# replace fields in 'x' with those in 'y', if given
## mergeRange <- function(x, y)
## {
##   if (length(genome(y)))
##     genome(x) <- genome(y)
##   if (length(y))
##     x[[1]] <- y[[1]]
##   if (length(names(y)))
##     names(x)[1] <- names(y)[1]
##   x
## }

## take 'y', filling in info from 'x' if necessary
mergeRange <- function(x, y)
{
  if (!length(genome(y))) {
    genome(y) <- genome(x)
    if (!length(names(y)) && length(y) == length(x)) {
      names(y) <- names(x)
      if (!length(unlist(y)) && length(x))
        y[[1]] <- x[[1]]
    }
  }
  if (!length(genome(y)))
    stop("Genome must be specified")
  if (length(unlist(y)) && !length(names(y)))
    stop("Chromosome must be specified if an interval is specified")
  y
}
