\name{feature-variables}
\alias{featChrom}
\alias{featStart}
\alias{featEnd}
\alias{featStrand}
\alias{featChrom,trackSet-method}
\alias{featStart,trackSet-method}
\alias{featEnd,trackSet-method}
\alias{featStrand,trackSet-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get track feature variables }
\description{
  Gets values for some of the most commonly used track feature variables
  (chrom, start, end, strand).
}
\usage{
featChrom(object)
featStart(object)
featEnd(object)
featStrand(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The \code{\linkS4class{trackSet}} holding the features. }
}
\value{
  The value of the variable (factor for chrom/strand, numeric for start/end).
}
\author{ Michael Lawrence }
\seealso{ \code{\link{trackData}} for getting all of these variables
  (and more) as a \code{data.frame}.}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
