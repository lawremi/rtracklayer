\name{browseGenome}
\alias{browseGenome}
\alias{browseGenome,ANY-method}
\title{ Browse a genome }
\description{
  A generic function for launching a genome browser.
}
\usage{
browseGenome(tracks = trackSets(), browser = "ucsc",
             segment = genomeSegment(tracks), view = TRUE,
             trackParams = list(), viewParams = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracks}{ A list of \code{\linkS4class{trackSet}} instances,
    e.g. a \code{\linkS4class{trackSets}} instance.}
  \item{browser}{ The name of the genome browser. }
  \item{segment}{ The \code{\linkS4class{genomeSegment}} to display in
    the initial view. }
  \item{view}{ Whether to open a view. }
  \item{trackParams} { Named list of parameters to pass to
    \code{\link{layTrack}}. }
  \item{viewParams} { Named list of parameters to pass to
    \code{\link{browserView}}. }
  \item{\dots}{ Arguments corresponding to slots in
    \code{\linkS4class{genomeSegment}} that override those in \code{segment}. }
}
\value{
  Returns a \code{\linkS4class{browserSession}}.
}
\author{ Michael Lawrence }
\seealso{ \code{\linkS4class{browserSession}} and
  \code{\linkS4class{browserView}}, the two main classes for interfacing
  with genome browsers.}
\examples{
  \dontrun{
  ## open UCSC genome browser:
  browseGenome()
  ## to view a specific segment:
  segment <- genomeSegment("hg18", "chr22", 20000, 50000)
  browseGenome(segment = segment)
  ## a slightly larger segment:
  browseGenome(segment = segment, end = 75000)
  ## with a track:
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  browseGenome(trackSets(track))
  }
}
\keyword{interface}
