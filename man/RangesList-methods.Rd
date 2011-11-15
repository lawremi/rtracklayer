\name{RangesList-methods}
\docType{methods}
\alias{RangesList-methods}

% accessors
\alias{chrom,RangesList-method}
\alias{chrom<-,RangesList-method}
\alias{seqinfo,RangesList-method}
\alias{seqinfo<-,RangesList-method}
\alias{seqnames,RangesList-method}

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
      \code{seqnames(x)}: Gets the sequence names for the ranges in
      \code{x}, as an \code{RleList} of length equal to the length
      of \code{x}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}: Gets or sets the
      sequence information as a
      \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} object. If this
      has not been set explicitly, it tries to come up with a reasonable
      default. First, it assumes the \code{\link[IRanges]{universe}} on
      \code{x} is a genome identifier and attempts to look up the
      corresponding metadata an installed BSgenome package or, if that
      fails, via UCSC. If that fails, it uses \code{names(x)} as the
      \code{seqnames} and \code{end(range(x))} as the
      \code{seqlengths}.
    }
  }
}

\author{ Michael Lawrence }

\keyword{classes}
\keyword{methods}
