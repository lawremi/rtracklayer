\name{IntegerRangesList-methods}
\docType{methods}
\alias{IntegerRangesList-methods}

% accessors
\alias{chrom,IntegerRangesList-method}
\alias{chrom<-,IntegerRangesList-method}

\title{Ranges on a Genome}

\description{
  Genomic coordinates are often specified in terms of a genome identifier,
  chromosome name, start position and end
  position. The \code{rtracklayer} package adds convenience methods to
  \code{IntegerRangesList} for the manipulation of genomic ranges.
  The spaces (or names) of \code{IntegerRangesList} are the chromosome names.
  The \code{universe} slot indicates the genome, usually as given by UCSC
  (e.g. \dQuote{hg18}).
}

\section{Accessors}{
  In the code snippets below,
  \code{x} is a \code{IntegerRangesList} object.

  \describe{
    \item{}{
      \code{chrom(x), chrom(x) <- value}: Gets or sets the chromosome
      names for \code{x}. This is an alias for
      \code{names(x)}.
    }
  }
}

\author{ Michael Lawrence }

\keyword{classes}
\keyword{methods}
