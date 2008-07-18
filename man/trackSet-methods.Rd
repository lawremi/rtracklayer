\name{trackSet-methods}
\docType{methods}
\alias{trackSet}
\alias{trackSet-methods}
\alias{trackSet,ucscSession-method}
\alias{trackSet,argoSession-method}
\title{Get a trackSet}
\description{
  Methods for getting a \code{\linkS4class{trackSet}}.
}
\usage{trackSet(object, ...)}
\arguments{
  \item{object}{An object from which to obtain a \code{trackSet}.}
  \item{\dots}{Arguments passed on to methods.}
}
\section{Methods}{
  \describe{
    The following methods are defined by \pkg{rtracklayer}.
    
    \item{object = "ucscSession"}{
      \code{trackSet(object, name, segment = genomeSegment(object), format =
        "bed", table = NULL)}:
      Get a track from the session by its \code{name} in format
      \code{format}. Some built-in tracks have multiple series, each
      identified by a \code{table} name.
    }
    
    \item{object = "argoSession"}{
      \code{trackSet(object, name, segment = genomeSegment(object))}:
      Get a track from the session by its \code{name}.
    }
  }
}
\seealso{
  \code{\link{layTrack}} for loading a \code{trackSet} into a session.
}
\examples{
\dontrun{
session <- browseGenome()
track <- import(system.file("tests", "v2.gff", package = "rtracklayer"))
layTrack(session, track, "My Track")
track <- trackSet(session, name = "My Track")
}
}
\keyword{methods}
