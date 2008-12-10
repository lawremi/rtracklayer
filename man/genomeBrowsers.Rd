\name{genomeBrowsers}
\alias{genomeBrowsers}
\title{Get available genome browsers}
\description{
  Gets the identifiers of the loaded genome browser drivers.
}
\usage{
genomeBrowsers(where = topenv(parent.frame()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{where}{ The environment in which to search for drivers. }
}
\details{
  This searches the specified environment for classes that extend
  \code{\linkS4class{BrowserSession}}. The prefix of the class name,
  e.g. "ucsc" in "UCSCSession", is returned for each driver. 
}
\value{
  A character vector of driver identifiers.
}
\author{ Michael Lawrence }
\seealso{ \code{\link{browseGenome}} and \code{\link{browserSession}}
  that create \code{browserSession} implementations given an identifier
  returned from this function. }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}
