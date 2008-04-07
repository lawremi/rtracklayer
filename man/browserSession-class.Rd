\name{browserSession-class}
\docType{class}
\alias{browserSession-class}
\alias{show,browserSession-method}

\title{Class "browserSession"}
\description{An object representing a genome browser session. Each
  session corresponds to a set of loaded \code{\linkS4class{trackSet}} instances
  and a set of \code{\linkS4class{browserView}} instances. Note that
  this is a virtual class; a concrete implementation is provided by each
  backend driver.}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it. See \code{\link{browserSession}} for obtaining
  an instance of an implementation for a particular genome browser.}
\section{Methods}{
  This specifies the API implemented by each browser backend. Note that
  a backend is not guaranteed to support all operations.
  \describe{
    \item{\code{\link{browserView}(object, segment =
        genomeSegment(object), track = tracks(object), ...)}}{
      Constructs a \code{\linkS4class{browserView}} of \code{segment}
      for this session.
    }
    \item{\code{\link{browserViews}(object, ...)}}{
      Gets the \code{\linkS4class{browserView}} instances belonging to
      this session.}
    \item{\code{\link{activeView}(object, ...)}}{
      Returns the \code{\linkS4class{browserView}} that is currently
      active in the session.}
    \item{\code{\link{genomeSegment}(object, ...)}}{
      Gets the \code{\linkS4class{genomeSegment}}
      representing the segment of the genome currently displayed by the
      browser (i.e. the segment shown by the active view) or a default
      value (possibly \code{NULL}) if no views exist.
    }
    \item{\code{\link{genomeSequence}(object, segment =
        genomeSegment(object), ...)}}{
      gets a genomic sequence of \code{segment} from this session.}
    \item{\code{\link{laySequence}(object, ...)}}{Loads a sequence into the
      session.}
    \item{\code{\link{layTrack}(object, track, name = deparse(substitute(track)),
        view = TRUE, ...)}}{
      Loads one or more \code{track}s into the session and optionally open a
      \code{view} of the track.
    }
    \item{\code{\link{trackSet}(object, segment = genomeSegment(object),
        name, ...)}}{
      Gets features in \code{segment} in the track stored as \code{name}
      as a \code{\linkS4class{trackSet}}.
    }
    \item{\code{\link{tracks}(object, ...)}}{ Gets the names of the
      tracks stored in this session. }
    \item{\code{\link{close}(con, ...)}}{ Close this session. }
    \item{\code{\link{show}(object, ...)}}{ Output a textual description
      of this session.}
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{browserSession}} for obtaining implementations of this
  class for a particular genome browser.
}
\keyword{classes}