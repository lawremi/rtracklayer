\name{GRangesForUCSCGenome}
\alias{GRangesForUCSCGenome}
\alias{GRangesForBSGenome}
\alias{SeqinfoForUCSCGenome}
\alias{SeqinfoForBSGenome}

\title{
  GRanges for a Genome
}
\description{
  These functions assist in the creation of
  \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} or
  \code{\link[GenomicRanges:GRanges-class]{GRanges}} for a genome.
}
\usage{
GRangesForUCSCGenome(genome, chrom = NULL, ranges = NULL, ...)
GRangesForBSGenome(genome, chrom = NULL, ranges = NULL, ...)

SeqinfoForUCSCGenome(genome)
SeqinfoForBSGenome(genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{
    A string identifying a genome, usually one assigned by UCSC, like "hg19".
  }
  \item{chrom}{
    A character vector of chromosome names, or \code{NULL}.
  }
  \item{ranges}{
    A \code{\link[IRanges:Ranges-class]{Ranges}} object with the intervals.
  }
  \item{\dots}{
    Additional arguments to pass to
    the \code{\link[GenomicRanges:GRanges]{GRanges}} constructor.
  }
}
\details{
  The genome ID is stored in the metadata of the ranges and is
  retrievable via the \code{\link{genome}} function. The sequence
  lengths are also properly initialized for the genome. This mitigates
  the possibility of accidentally storing intervals for the wrong
  genome.

  \code{GRangesForUCSCGenome} obtains sequence information from the UCSC
  website, while \code{GRangesForBSGenome} looks for it in an
  installed \code{BSGenome} package. Using the latter is more efficient
  in the long-run, but requires downloading and installing a potentially
  large genome package, or creating one from scratch if it does not yet
  exist for the genome of interest.
}
\value{
  For the \code{GRangesFor*} functions, a \code{GRanges} object, with the
  appropriate \code{\link[GenomicRanges:seqlengths]{seqlengths}} and
  \code{\link{genome}} ID.

  The \code{SeqinfoFor*} functions return a \code{Seqinfo} for the
  indicated genome.
}
\author{
  Michael Lawrence
}
