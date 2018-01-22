\name{GenomicData}

% constructor
\alias{GenomicData}

% accessors
\alias{chrom}
\alias{chrom,GRanges-method}
\alias{chrom<-}
\alias{chrom<-,GRanges-method}
\alias{score,ANY-method}

\title{Data on a Genome}

\description{The \code{rtracklayer} package adds
convenience methods on top of \code{GenomicRanges} and \code{IntegerRangesList}
to manipulate data on genomic ranges.}

\section{Accessors}{
  In the code snippets below,
  \code{x} is a \code{GenomicRanges} or \code{IntegerRangesList} object.

  \describe{
    \item{}{
      \code{chrom(x), chrom(x) <- value}: Gets or
      sets the chromosome names for \code{x}. The length of \code{value}
      should equal the length of \code{x}.
    }
    \item{}{
      \code{score(x)}: Gets the \dQuote{score} column from the element
      metadata of a \code{GenomicRanges} or \code{GRangesList}. Many
      track formats have a score column, so this is often used during
      export. The \code{ANY} fallback for this method simply returns
      \code{NULL}.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicData(ranges, ..., strand = NULL, chrom = NULL,
        genome = NULL)}: Constructs a \code{GRanges} instance with
      the given \code{ranges} and variables in \code{...} (see the
      \code{\link[GenomicRanges]{GRanges}} constructor).

      If non-\code{NULL}, the \code{strand} argument specifies the strand
      of each range. It should be a character vector or factor of length
      equal to that of \code{ranges}. All values should be either \code{-},
      \code{+}, or \code{*}. To get the levels for \code{strand}, call
      \code{levels(strand())}.

      \code{chrom} argument is analogous to \code{seqnames} in the
      \code{GRanges} constructor.

      The \code{genome} argument should be a scalar string. See the examples.
    }
  }
}

\author{ Michael Lawrence and Patrick Aboyoun }

\examples{
  range1 <- IRanges(c(1,2,3), c(5,2,8))

  ## with some data ##
  filter <- c(1L, 0L, 1L)
  score <- c(10L, 2L, NA)
  strand <- factor(c("+", NA, "-"), levels = levels(strand()))
  ## GRanges instance
  gr <- GenomicData(range1, score, chrom = "chr1", genome = "hg18")
  mcols(gr)[["score"]]
  strand(gr) ## all '*'
  gr <- GenomicData(range1, score, filt = filter, strand = strand,
                    chrom = "chr1")
  mcols(gr)[["filt"]]
  strand(gr) ## equal to 'strand'
  ## coercion from data.frame ##
  df <- as.data.frame(gr)
}
\keyword{classes}
\keyword{methods}
