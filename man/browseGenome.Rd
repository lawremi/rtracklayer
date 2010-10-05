\name{browseGenome}
\alias{browseGenome}
\alias{browseGenome,RangedDataORRangedDataList-method}
\alias{browseGenome,missing-method}
\title{ Browse a genome }
\description{
  A generic function for launching a genome browser.
}
\usage{
browseGenome(object, ...)
\S4method{browseGenome}{RangedDataORRangedDataList}(object,
  browser = "UCSC", range = base::range(object),
  view = TRUE, trackParams = list(), viewParams = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A list of
    \code{\link[IRanges:RangedData-class]{RangedData}} instances, 
    e.g. a \code{\link[IRanges:RangedDataList-class]{RangedDataList}} instance.}
  \item{browser}{ The name of the genome browser. }
  \item{range}{ A genome identifier or a
    \code{\link[GenomicRanges:GRanges-class]{GRanges}} to display in the
    initial view.
    Support for \code{\link[IRanges:RangesList-class]{RangesList}} is also
    available, but deprecated. }
  \item{view}{ Whether to open a view. }
  \item{trackParams}{ Named list of parameters to pass to
    \code{\link{track<-}}. }
  \item{viewParams}{ Named list of parameters to pass to
    \code{\link{browserView}}. }
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
  range <- GenomicRanges("hg18", "chr22", 20000, 50000)
  browseGenome(range = range)
  ## a slightly larger range:
  browseGenome(range = range, end = 75000)
  ## with a track:
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  browseGenome(RangedDataList(track))
  }
}
\keyword{interface}
