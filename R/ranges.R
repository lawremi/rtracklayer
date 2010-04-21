### =========================================================================
### Genome-oriented methods for RangedData classes
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("genome", function(x, ...) standardGeneric("genome"))
setMethod("genome", "RangedData", function(x) universe(x))

setGeneric("genome<-", function(x, value) standardGeneric("genome<-"))
setReplaceMethod("genome", "RangedData",
                 function(x, value) {
                   genome(x@ranges) <- value
                   x
                 })

setGeneric("chrom", function(x, ...) standardGeneric("chrom"))
setMethod("chrom", "RangedData", function(x) {
  chrom(ranges(x))
})

setGeneric("chrom<-", function(x, ..., value) standardGeneric("chrom<-"))
setReplaceMethod("chrom", "RangedData", function(x, value) {
  chrom(ranges(x)) <- value
  x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicData <- function(ranges, ..., strand = NULL, chrom = NULL, genome = NULL)
{
  if (is(ranges, "data.frame") || is(ranges, "DataTable")) {
    colnames(ranges)[match("chrom", colnames(ranges))] <- "space"
  }
  if (!is(ranges, "Ranges")) {
    return(RangedData(ranges)) # direct coercion
  }
  if (length(chrom) > length(ranges))
    stop("length of 'chrom' greater than length of 'ranges'")
  if (length(chrom) > 0 && (length(ranges) %% length(chrom) != 0))
    stop("length of 'ranges' not a multiple of 'chrom' length")
  if (!is.null(genome) && (length(genome) != 1 || !is.character(genome)))
    stop("'genome' must be a single string")
  if (!is.null(strand)) {
    if (!all(strand[!is.na(strand)] %in% levels(strand())))
      stop("strand values should be 'NA', '-', '+' or '*'")
    strand <- factor(strand, levels(strand()))
    RangedData(ranges, ..., strand = strand, space = chrom, universe = genome)
  } else RangedData(ranges, ..., space = chrom, universe = genome)
}

### =========================================================================
### Genome-oriented methods for RangesList classes
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

GenomicRanges <- function(start = integer(), end = integer(), chrom = NULL,
                          genome = NULL)
{
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
    stop("'genome' does not correspond to an installed BSgenome package")
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
