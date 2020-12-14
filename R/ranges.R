### =========================================================================
### Genome-oriented methods for GRanges and IntegerRangesList objects
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GRanges convenience constructor
###

GenomicData <- function(ranges, ..., strand = NULL, chrom = NULL, genome = NA,
                        seqinfo = NULL, which = NULL, metadata = list())
{
  if (is.null(genome))
    genome <- NA
  if (!isSingleStringOrNA(genome))
    stop("'genome' must be a single string, or NULL, or NA")
  if (!is.null(seqinfo) && !all(genome == genome(seqinfo), na.rm=TRUE)) {
      stop("'genome' ", genome, "' does not match that in 'seqinfo'")
  }
  if (length(seqinfo) == 0L)
    seqinfo <- seqinfoForGenome(genome)
  if (length(chrom) > length(ranges))
    stop("length of 'chrom' greater than length of 'ranges'")
  if (length(chrom) > 0 && (length(ranges) %% length(chrom) != 0))
    stop("length of 'ranges' not a multiple of 'chrom' length")
  normStrand <- function(strand) {
    strand <- as.character(strand)
    strand[is.na(strand)] <- "*"
    strand(strand)
  }
  if (!(is.null(strand) || is(strand, "Rle")))
    strand <- normStrand(strand)
  if (is.null(chrom))
    chrom <- droplevels(Rle(factor("1"), length(ranges)))
  dots <- list(...)
  if (length(dots) == 1) {
    dots <- dots[[1L]]
    if ((is(dots, "data.frame") || is(dots, "DataFrame")) &&
        !is.null(dots[["strand"]])) {
      strand <- dots[["strand"]]
      dots[["strand"]] <- NULL
      return(GenomicData(ranges = ranges, dots, strand = strand,
                         chrom = chrom, genome = genome,
                         seqinfo = seqinfo, which = which,
                         metadata = metadata))
    }
  }
  if (is.null(strand))
    strand <- Rle("*", length(ranges))
  df <- DataFrame(...)
  invalidNames <- names(df) %in% GenomicRanges:::INVALID.GR.COLNAMES
  names(df)[invalidNames] <- paste0(".", names(df)[invalidNames])
  gd <- GRanges(seqnames = chrom, ranges = ranges, strand = strand, df,
                seqinfo = seqinfo)
  if (!is.na(genome))
    genome(gd) <- genome
  if (!is.null(which)) {
    if (is(which, "IntegerRangesList"))
      which <- as(which, "GRanges")
    gd <- subsetByOverlaps(gd, which)
  }
  metadata(gd) <- metadata
  gd
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Automatic seqinfo lookup
###

seqinfoForGenome <- function(genome, method = c("auto", "BSgenome", "UCSC")) {
  method <- match.arg(method)
  if (is.na(genome))
    return(NULL)
  if (method == "auto" || method == "BSgenome")
    sl <- SeqinfoForBSGenome(genome)
  if (method == "UCSC" || (method == "auto" && is.null(sl)))
    sl <- SeqinfoForUCSCGenome(genome)
  sl
}

BSGenomeForID <- function(genome) {
  if (!suppressWarnings(requireNamespace("BSgenome", quietly=TRUE)))
    return(NULL)
  bsgenome <- try(BSgenome::getBSgenome(genome), silent=TRUE)
  if (inherits(bsgenome, "try-error"))
    return(NULL)
  bsgenome
}

SeqinfoForBSGenome <- function(genome) {
  bsgenome <- BSGenomeForID(genome)
  if (!is.null(bsgenome))
    seqinfo(bsgenome)
  else NULL
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Convience constructor for GRanges over an entire genome
###

GRangesForGenome <- function(genome, chrom = NULL, ranges = NULL, ...,
                             method = c("auto", "BSgenome", "UCSC"),
                             seqinfo = NULL)
{
  if (missing(genome) || !isSingleString(genome))
    stop("'genome' must be a single string identifying a genome")
  if (is.null(seqinfo))
    seqinfo <- seqinfoForGenome(genome, match.arg(method))
  if (is.null(seqinfo))
    stop("Failed to obtain information for genome '", genome, "'")
  if (!is.null(ranges) && !is(ranges, "IntegerRanges"))
    stop("'ranges' must be NULL or an IntegerRanges object")
  if (is.null(chrom))
    chrom <- seqnames(seqinfo)
  else {
    badChrom <- setdiff(chrom, seqnames(seqinfo))
    if (length(badChrom))
      stop("Chromosome(s) ", paste(badChrom, collapse = ", "),
           "are invalid for: ", genome)
  }
  if (is.null(ranges))
    ranges <- IRanges(1L, seqlengths(seqinfo)[chrom])
  gr <- GRanges(chrom, ranges, seqlengths = seqlengths(seqinfo), ...)
  seqinfo(gr) <- seqinfo
  gr
}

GRangesForBSGenome <- function(genome, chrom = NULL, ranges = NULL, ...)
{
  GRangesForGenome(genome, chrom = chrom, ranges = ranges, method = "BSgenome",
                   seqinfo = NULL, ...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### score: for internal convenience
###

setMethod("score", "ANY", function(x) NULL)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chrom(): Returns chromosome name vector of length 'length(x)'.
###          Not to be confused with 'seqnames', which returns a List for
###          IntegerRangesList.
###          More or less a pre-GenomicRanges relic.
###

setGeneric("chrom", function(x, ...) standardGeneric("chrom"))
setMethod("chrom", "GRanges", function(x) seqnames(x))
setMethod("chrom", "IntegerRangesList", function(x) {
  names(x)
})

setGeneric("chrom<-", function(x, ..., value) standardGeneric("chrom<-"))
setReplaceMethod("chrom", "GRanges", function(x, value) {
  seqnames(x) <- value
  x
})

### =========================================================================
### Genome-oriented conveniences for RangedSelection classes
### -------------------------------------------------------------------------

GenomicSelection <- function(genome, chrom = NULL, colnames = character(0))
{
  if (missing(genome) || !isSingleString(genome))
    stop("'genome' must be a single string identifying a genome")
  si <- seqinfoForGenome(genome)
  if (is.null(si))
    stop("Failed to obtain information for genome '", genome, "'")
  lens <- seqlengths(si)
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
### Utilities
###

## Individual sessions/views only support a single genome at a time,
## whereas range data structures can have multiple genomes. We try to
## rectify this here.
singleGenome <- function(x) {
  x1 <- unname(x[1])
  if (!all(is.na(x)) && (any(is.na(x)) || any(x != x1)))
      stop("Multiple genomes encountered; only one supported")
  x1
}

## normalize 'range', using 'seqinfo' for default genome
normGenomeRange <- function(range, seqinfo, max.length = 1L) {
  ## the user can specify a portion of the genome in several ways:
  ## - String identifying a genome
  ## - IntegerRangesList
  ## - GRanges, the preferred way, possibly from GRangesForUCSCGenome()
  ## - We do not allow IntegerRanges, since it does not make sense to have
  ##   one range over many chromosomes
  if (is.character(range)) {
    range <- singleGenome(range)
    seqinfo <- Seqinfo(genome = range)
    return(GRangesForGenome(range, seqinfo = seqinfo))
  }
  genome <- genome(seqinfo)[1]
  if (is(range, "Seqinfo"))
    range <- as(range, "GRanges")
  if (!is(range, "IntegerRangesList") && !is(range, "GenomicRanges"))
    stop("'range' should be a genome string, IntegerRangesList, GRanges or Seqinfo")

  if (length(seqinfo(range)) == 0L) {
    ## hack: need to avoid calling seqlengths(session) here, so use 'foo'
    seqinfo(range) <- Seqinfo("foo", genome = genome)
  } else {
    rangeGenome <- singleGenome(genome(range))
    if (!is.na(rangeGenome) && rangeGenome != genome) {
      genome <- rangeGenome
    }
    seqinfo <- si <- Seqinfo(genome = genome)
    seqinfo(range, new2old = match(seqlevels(si), seqlevels(range))) <-
      merge(si, seqinfo(range))
  }
  if (is(range, "IntegerRangesList")) {
    range <- GRangesForGenome(singleGenome(genome(range)), names(range),
                              unlist(range), seqinfo = seqinfo)
  } else if (is(range, "GenomicRanges")) {
    strand(range) <- "*"
    mcols(range) <- NULL
  }
  if (length(range) > max.length) {
    warning("number of ranges (", length(range),
            ") exceeds limit of ", max.length)
  }
  if (max.length == 1L) {
    if (length(unique(seqnames(range))) > max.length)
      stop("'range' must contain ranges on a single chromosome")
    range(range)
  } else {
    head(range, max.length)
  }
}

spansGenome <- function(x) {
  strand(x) <- "*"
  x <- reduce(x)
  w <- structure(width(x), names = as.character(seqnames(x)))
  identical(w[names(seqlengths(x))], seqlengths(x))
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
