\name{BigWigSelection-class}
\docType{class}
\alias{BigWigSelection-class}

% constructor
\alias{BigWigSelection}

\alias{coerce,RangesList,BigWigSelection-method}
\alias{coerce,GenomicRanges,BigWigSelection-method}

\title{Selection of ranges and columns}

\description{A \code{BigWigSelection} represents a query against a
  BigWig file, see \code{\link{import.bw}}. It is simply
  a \link[IRanges]{RangedSelection} that requires its \code{colnames}
  parameter to be "score", if non-empty, as that is the only column
  supported by BigWig.}

\section{Constructor}{
  \describe{
    \item{}{\code{BigWigSelection(ranges = GRanges(), colnames =
        "score")}: Constructs a \code{BigWigSelection} with the given
        \code{ranges} and \code{colnames}. \code{ranges} can be either
        something coercible to a \code{\linkS4class{RangesList}}, a
        \code{character} identifying a genome (see
        \code{\link{GenomicSelection}}), or a
        \code{\linkS4class{BigWigFile}}, in which case the ranges are
        derived from the bounds of its sequences.
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{\code{as(from, "BigWigSelection")}: Coerces \code{from} to a
      \code{BigWigSelection} object. Typically, \code{from} is a
      \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
      a \code{\link[IRanges:RangesList-class]{RangesList}}, the ranges of
      which become the ranges in the
      new \code{BigWigSelection}.
    }
  }
}

\author{ Michael Lawrence }

\examples{
  rl <- IRanges::RangesList(chr1 = IRanges::IRanges(c(1, 5), c(3, 6)))

  BigWigSelection(rl)
  as(rl, "BigWigSelection") # same as above

  # do not select the 'score' column
  BigWigSelection(rl, character())
}

\keyword{methods}
\keyword{classes}
