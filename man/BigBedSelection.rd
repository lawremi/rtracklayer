\name{BigBedSelection-class}
\docType{class}
\alias{BigBedSelection-class}

% constructor
\alias{BigBedSelection}

\alias{coerce,IntegerRangesList,BigBedSelection-method}
\alias{coerce,GenomicRanges,BigBedSelection-method}

\title{Selection of ranges and columns}

\description{A \code{BigBedSelection} represents a query against a
  BigBed file, see \code{\link{import.bb}}. It is simply
  a \link[IRanges]{RangedSelection} with \code{colnames}
  parameter.\code{colnames} should be a character vector of column names.
  Default columns are \code{"name", "score", "thick", "itemRgb"}
  and \code{"blocks"}, if non-empty, as that is the only column supported
  by BigBed.}

\section{Constructor}{
  \describe{
    \item{}{\code{BigBedSelection(ranges = GRanges(), colnames =
        "score")}: Constructs a \code{BigBedSelection} with the given
        \code{ranges} and \code{colnames}.
        a \code{character} identifying a genome (see
        \code{\link{GenomicSelection}}), or a
        \code{\linkS4class{BigBedFile}}, in which case the ranges are
        derived from the bounds of its sequences.
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{\code{as(from, "BigBedSelection")}: Coerces \code{from} to a
      \code{BigBedSelection} object. Typically, \code{from} is a
      \code{\link[GenomicRanges]{GRanges}} or
      a \code{\link[IRanges]{IntegerRangesList}}, the ranges of
      which become the ranges in the
      new \code{BigBedSelection}.
    }
  }
}

\author{ Michael Lawrence }

\examples{
  rl <- IRangesList(chr1 = IRanges::IRanges(c(1, 5), c(3, 6)))

  BigBedSelection(rl)
  as(rl, "BigBedSelection") # same as above

  # do not select any column
  BigBedSelection(rl, character())
}

\keyword{methods}
\keyword{classes}
