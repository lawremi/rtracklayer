### =========================================================================
### Genome-oriented methods for GRanges/RangedData classes
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
              metadata(x)$universe # 'universe' for compat with RangedData
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
                   metadata(x)$universe <- value
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

## zooming (symmetrically scales the width)
setMethod("Ops", c("GRanges", "numeric"),
          function(e1, e2)
          {
            if (IRanges:::anyMissing(e2))
              stop("NA not allowed as zoom factor")
            if ((length(e1) < length(e2) && length(e1)) ||
                (length(e1) && !length(e2)) ||
                (length(e1) %% length(e2) != 0))
              stop("zoom factor length not a multiple of number of ranges")
            if (.Generic == "*") {
              e2 <- ifelse(e2 < 0, abs(1/e2), e2)
              resize(e1, width(e1) / e2, fix = "center")
            } else {
              if (.Generic == "-") {
                e2 <- -e2
                .Generic <- "+"
              }
              if (.Generic == "+") {
                if (any(-e2*2 > width(e1)))
                  stop("adjustment would result in ranges with negative widths")
                resize(e1, width(e1) + e2*2, fix = "center")
              }
            }
          }
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

seqinfoForGenome <- function(genome, method = c("auto", "BSgenome", "UCSC")) {
  method <- match.arg(method)
  if (method == "auto" || method == "BSgenome")
    sl <- seqinfoForBSGenome(genome)
  if (method == "UCSC" || (method == "auto" && is.null(sl)))
    sl <- seqinfoForUCSCGenome(genome)
  sl
}

seqinfoForBSGenome <- function(genome) {
  bsgenome <- BSGenomeForID(genome)
  if (!is.null(bsgenome))
    seqinfo(bsgenome)
  else NULL
}

GRangesForGenome <- function(genome, chrom = NULL, ranges = NULL,
                             method = c("auto", "BSgenome", "UCSC"), ...)
{
  if (missing(genome) || !IRanges:::isSingleString(genome))
    stop("'genome' must be a single string identifying a genome")
  si <- seqinfoForGenome(genome, match.arg(method))
  if (is.null(si))
    stop("Failed to obtain information for genome '", genome, "'")
  if (!is.null(ranges) && !is(ranges, "Ranges"))
    stop("'ranges' must be NULL or a Ranges object")
  if (is.null(chrom))
    chrom <- seqnames(si)
  else {
    badChrom <- setdiff(chrom, seqnames(si))
    if (length(badChrom))
      stop("Chromosome(s) ", paste(badChrom, collapse = ", "),
           "are invalid for: ", genome)
  }
  if (is.null(ranges))
    ranges <- IRanges(1L, seqlengths(si)[chrom])
  gr <- GRanges(chrom, ranges, seqlengths = seqlengths(si), ...)
  seqinfo(gr) <- si
  genome(gr) <- genome
  gr
}

GRangesForBSGenome <- function(genome, chrom = NULL, ranges = NULL, ...)
{
  GRangesForGenome(genome, chrom = chrom, ranges = ranges, method = "BSgenome",
                   ...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

## normalize 'range', using 'session' for default genome
## note that is singular, i.e., only one interval should come out of this
normGenomeRange <- function(range, session) {
  ## the user can specify a portion of the genome in several ways:
  ## - String identifying a genome
  ## - RangesList, possibly constructed with deprecated GenomicRanges()
  ## - GRanges, the preferred way, possibly from GRangesForUCSCGenome()
  ## - We do not allow Ranges, since it does not make sense to have one range
  ##   over many chromosomes
  if (is.character(range)) {
    return(GRangesForUCSCGenome(range))
  }
  genome <- genome(session)  
  if (is.null(genome(range)))
    genome(range) <- genome
  else if (genome(range) != genome) {
    genome(session) <- genome(range)
    on.exit(genome(session) <- genome)
  }
  if (is(range, "RangesList")) {
    if (length(range) > 1)
      stop("'range' must contain ranges on a single chromosome")
    chrom <- names(range)
    genome <- genome(range)
    ranges <- NULL
    if (!is.null(chrom)) {
      flatRange <- unlist(range)
      if (length(flatRange))
        range <- range(flatRange)
    }
    GRangesForUCSCGenome(genome, chrom, range)
  } else if (is(range, "GRanges")) {
    if (length(unique(seqnames(range))) != 1L)
      stop("'range' must contain ranges on a single chromosome")
    strand(range) <- "*"
    range <- range(range)
    seqlens <- seqlengths(session)
    chr <- as.character(seqnames(range))
    if (!chr %in% names(seqlens))
      stop("'range' has invalid chromosome: ", chr)
    if (start(range) < 1L || end(range) > seqlens[chr])
      stop("'range' is out of bounds for ", genome(range), ":", chr)
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

BSGenomeForID <- function(genome) {
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
