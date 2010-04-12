\name{BigWigSelection-class}
\docType{class}
\alias{BigWigSelection-class}

% constructor
\alias{BigWigSelection}

\alias{coerce,RangesList,BigWigSelection-method}

\title{Selection of ranges and columns}

\description{A \code{BigWigSelection} represents a query against a
  BigWig file, see \code{\link{import.bw}}. It is simply
  a \link[IRanges]{RangedSelection} that requires its \code{colnames}
  parameter to be "score", if non-empty, as that is the only column
  supported by BigWig.}

\section{Constructor}{
  \describe{
    \item{}{\code{BigWigSelection(ranges = RangesList(),
        colnames = "score")}: Constructors a \code{BigWigSelection}
        with the given \code{ranges} and \code{colnames}.
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{\code{as(from, "BigWigSelection")}: Coerces \code{from} to a
      \code{BigWigSelection} object. Typically, \code{from} is a
      \code{\linkS4class{RangesList}}, the ranges of which become the
      ranges in the new \code{BigWigSelection}.
    }
  }
}

\author{ Michael Lawrence }

\examples{
  rl <- RangesList(chr1 = IRanges(c(1, 5), c(3, 6)))

  BigWigSelection(rl)
  as(rl, "BigWigSelection") # same as above

  # do not select the 'score' column
  BigWigSelection(rl, character())
}

\keyword{methods}
\keyword{classes}