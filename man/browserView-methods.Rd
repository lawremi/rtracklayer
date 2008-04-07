\name{browserView-methods}
\docType{methods}
\alias{browserView}
\alias{browserView-methods}
\alias{browserView,ucscSession-method}
\alias{browserView,argoSession-method}
\title{Getting browser views}
\description{
  Methods for creating and getting browser views.
}
\usage{
  browserView(object, segment = genomeSegment(object),
              track = tracks(object), ...)
}
\arguments{
  \item{object}{The object from which to get the views.}
  \item{segment}{The \code{\linkS4class{genomeSegment}} to display.}
  \item{track}{List of track names to make visible in the view.}
  \item{\dots}{Arguments to pass to methods}
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
  \describe{
    
    \item{object = "ucscSession"}{ \code{browserView(object, segment =
        genomeSegment(object), track = tracks(object), ...)}: Creates a
      \code{\linkS4class{browserView}} of \code{segment} with visible
      tracks named in \code{track}. Arguments in \code{...} correspond
      to slots in the \code{\linkS4class{genomeSegment}} class and override those
      in \code{segment}. }
    
    \item{object = "argoSession"}{ \code{browserView(object, segment =
        genomeSegment(object), track = tracks(object, segment, TRUE),
        ...)}: Creates a \code{\linkS4class{browserView}} of
      \code{segment} with visible tracks named in
      \code{track}. Parameters in \code{...} correspond to slots in
      the \code{\linkS4class{genomeSegment}} class and override those in
      \code{segment}. }
}}
\examples{
\dontrun{
  session <- browserSession()
  browserView(session, genomeSegment(start = 20000, end = 50000,
                                     segment = genomeSegment(session)))
  ## equivalent to above, but shorter
  browserView(session, start = 20000, end = 50000)
  ## only view "knownGene" track
  browserView(session, track = "knownGene")
}
}
\keyword{methods}
