\name{RangesList-methods}
\docType{methods}
\alias{RangesList-methods}

% constructor
\alias{GenomicRanges}

% accessors
\alias{chrom,RangesList-method}
\alias{chrom<-,RangesList-method}
\alias{genome,RangesList-method}
\alias{genome<-,RangesList-method}

\title{Ranges on a Genome}

\description{
  Genomic coordinates are often specified in terms of a genome identifier,
  chromosome name, start position and end
  position. \code{\link[IRanges:RangedData-class]{RangedData}} represents
  this with a \code{\link[IRanges:RangesList-class]{RangesList}}
  instance, and the \code{rtracklayer} package adds convenience methods to
  \code{RangesList} for the manipulation of genomic ranges. The spaces (or
  names) of \code{RangesList} are the chromosome names. The
  \code{universe} slot indicates the genome, usually as given by UCSC
  (e.g. \dQuote{hg18}).
}

\section{Accessors}{
  In the code snippets below,
  \code{x} is a \code{RangesList} object.

  \describe{
    \item{}{
      \code{chrom(x), chrom(x) <- value}: Gets or sets the chromosome
      names for \code{x}. This is an alias for
      \code{\link[IRanges:RangesList-class]{names}(x)}.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}: Gets or sets the
      genome (a single string or \code{NULL}) for
      the ranges in \code{x}; simple wrappers around
      \code{\link[IRanges:RangesList-class]{universe}} and
      \code{\link[IRanges:RangesList-class]{universe<-}}, respectively.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicRanges(start, end, chrom = NULL, genome = NULL)}:
      
      \strong{This function has been deprecated in favor of constructing
        a \code{\link[GenomicRanges:GRanges-class]{GRanges}} using
        \code{\link{GRangesForUCSCGenome}}.}
      
      Constructs a \code{RangesList} containing ranges specified by
      \code{start} and \code{end}, optionally split into elements based on
      \code{chrom}, a vector of chromosome identifiers (or \code{NULL} for no
      splitting). The \code{genome} argument should be a
      scalar string and is treated as the \code{RangesList}
      universe. See the examples.
    }
  }
}

\author{ Michael Lawrence }

\examples{
GenomicRanges(c(1,2,3), c(5,2,8))
GenomicRanges(c(1,2,3), c(5,2,8), c("chr1", "chr1", "chr2"))
GenomicRanges(c(1,2,3), c(5,2,8), genome = "hg18")
}
\keyword{classes}
\keyword{methods}
