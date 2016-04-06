\name{RangesList-methods}
\docType{methods}
\alias{RangesList-methods}

% accessors
\alias{chrom,RangesList-method}
\alias{chrom<-,RangesList-method}

\title{Ranges on a Genome}

\description{
  Genomic coordinates are often specified in terms of a genome identifier,
  chromosome name, start position and end
  position. The \code{rtracklayer} package adds convenience methods to
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
  }
}

\author{ Michael Lawrence }

\keyword{classes}
\keyword{methods}
