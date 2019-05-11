\name{browseGenome}
\alias{browseGenome}
\alias{browseGenome,GenomicRanges_OR_GenomicRangesList-method}
\alias{browseGenome,missing-method}
\title{ Browse a genome }
\description{
  A generic function for launching a genome browser.
}
\usage{
browseGenome(object, ...)
\S4method{browseGenome}{GenomicRanges_OR_GenomicRangesList}(object,
  browser = "UCSC", range = base::range(object),
  view = TRUE, trackParams = list(), viewParams = list(),
  name = "customTrack", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{\link[GenomicRanges]{GRanges}} object or a list
    of \code{\link[GenomicRanges]{GRanges}} objects (e.g. a
    \code{\link[GenomicRanges]{GenomicRangesList}} object). }
  \item{browser}{ The name of the genome browser. }
  \item{range}{ A genome identifier or a
    \code{\link[GenomicRanges]{GRanges}} or
  \code{\link[IRanges]{IntegerRangesList}} to display in the
    initial view. }
  \item{view}{ Whether to open a view. }
  \item{trackParams}{ Named list of parameters to pass to
    \code{\link{track<-}}. }
  \item{viewParams}{ Named list of parameters to pass to
    \code{\link{browserView}}. }
  \item{name}{ The name for the track. }
  \item{\dots}{ Arguments passed to \code{\link{browserSession}}. } 
}
\value{
  Returns a \code{\linkS4class{BrowserSession}}.
}
\author{ Michael Lawrence }
\seealso{ \code{\linkS4class{BrowserSession}} and
  \code{\linkS4class{BrowserView}}, the two main classes for interfacing
  with genome browsers.}
\examples{
  \dontrun{
  ## open UCSC genome browser:
  browseGenome()
  ## to view a specific range:
  range <- GRangesForUCSCGenome("hg18", "chr22", IRanges(20000, 50000))
  browseGenome(range = range)
  ## a slightly larger range:
  browseGenome(range = range, end = 75000)
  ## with a track:
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  browseGenome(GRangesList(track))
  }
}
\keyword{interface}
