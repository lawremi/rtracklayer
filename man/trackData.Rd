\name{trackData}
\alias{trackData}
\alias{trackData,trackSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get track data }
\description{
  Convert the track features and experimental measurements to a
  \code{data.frame}.
}
\usage{
trackData(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object with track information, e.g.
    \code{\linkS4class{trackSet}}. }
  \item{\dots}{ Arguments to pass to methods }
}
\value{
  A \code{data.frame} with a column for each feature variable and
  experimental sample. An additional column \code{seqMid} is added that
  holds the mid-point of each feature (helpful for plotting).
}
\author{ Michael Lawrence }
\seealso{ \code{\link{featChrom}}, \code{\link{featStart}},
  \code{\link{dataVals}}, etc for accessing individual elements of the
  track.}
\keyword{ manip }
