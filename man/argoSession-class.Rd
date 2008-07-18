\name{argoSession-class}
\docType{class}
\alias{argoSession-class}
\alias{close,argoSession-method}
\alias{initialize,argoSession-method}

\title{Class "argoSession"}
\description{An object representing a session in the Argo genome browser.}
\section{Objects from the Class}{
  Objects may be created by calls of the form
  \code{\link{browserSession}("argo", name, jar)} where \code{name} is
  the name of the session in Argo and \code{jar} is the path to the
  installed Argo jar file.
}
\section{Slots}{
  \describe{
    \item{\code{jref}:}{Object of class \code{"jobjRef"} holding the
      reference to the low-level Argo object. }
    \item{\code{seqcache}:}{Object of class \code{"environment"} that
      caches sequence objects after lookup from \acronym{UCSC}. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{browserSession}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{activeView}(object)}}{ Obtains the currently
      active \code{\linkS4class{browserView}} for this session. }
    \item{\code{\link{activeView<-}(object)}}{ Sets the currently active
      view for this session.}
    \item{\code{\link{browserView}(object, segment =
        genomeSegment(object),track = tracks(object, segment, TRUE),
        ...)}}{
      Creates a \code{\linkS4class{browserView}} of
      \code{segment} with visible tracks named in
      \code{track}. Parameters in \code{...} correspond to slots in
      the \code{\linkS4class{genomeSegment}} class and override those in
      \code{segment}.}
    \item{\code{\link{browserViews}(object)}}{ Gets a list of the
      \code{\linkS4class{browserView}} instances for this session. }
    \item{\code{\link{close}(con)}}{ Closes this session. }
    \item{\code{\link{genomeSequence}(object, segment)}}{ Gets the
      genome sequence in \code{segment}. }
    \item{\code{\link{laySequence}(object, sequence, name, label =
        name)}}{ Stores \code{sequence} under \code{name}, labeled as
      \code{label} in the user interface. }
    \item{\code{\link{layTrack}(object, track, name =
        deparse(substitute(track)), view = TRUE)}}{ Stores \code{track}
      under \code{name}, opening a view of the track if \code{view} is
      \code{TRUE}.}
    \item{\code{\link{trackSet}(object, name,
        segment = genomeSegment(object))}}{
      Obtains a \code{\linkS4class{trackSet}} named \code{name} from
      \code{segment}.}
    \item{\code{\link{tracks}(object, segment = NULL, visible =
        FALSE)}}{ Gets the names of tracks in \code{segment}. If
      \code{visible} is \code{TRUE}, only the visible tracks are
      returned. Note that in Argo tracks are visible on a
      global, not per-session nor per-view, basis. }
  }
}
\references{The Argo genome browser:
  \url{http://www.broad.mit.edu/annotation/argo/}}

\author{ Michael Lawrence }
\note{ In order to use this backend, the Argo jar file must be manually
  downloaded from the \acronym{URL} in the references. The call to
  \code{\link{browserSession}} must then pass the path to the jar file
  as the \code{jar} parameter. By default, the path is the \file{java}
  directory of the \pkg{rtracklayer} package. }

\seealso{
  \code{\link{browserSession}} for creating instances of this class.
}
\keyword{classes}
