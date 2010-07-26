\name{BrowserSession-class}
\docType{class}
\alias{BrowserSession-class}
\alias{show,BrowserSession-method}
\alias{[[,BrowserSession-method}
\alias{[[<-,BrowserSession-method}
\alias{$,BrowserSession-method}
\alias{$<-,BrowserSession-method}
\alias{range,BrowserSession-method}
\alias{genome,BrowserSession-method}
\alias{genome<-,BrowserSession-method}


\title{Class "BrowserSession"}
\description{An object representing a genome browser session. Each
  session corresponds to a set of loaded tracks
  and a set of \code{\linkS4class{BrowserView}} instances. Note that
  this is a virtual class; a concrete implementation is provided by each
  backend driver.}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it. See \code{\link{browserSession}} for obtaining
  an instance of an implementation for a particular genome browser.}
\section{Methods}{
  This specifies the API implemented by each browser backend. Note that
  a backend is not guaranteed to support all operations, and that each
  backend often has additional parameters for each of the methods. See
  the backend-specific documentation for more details. The only built-in
  backend is \code{\linkS4class{UCSCSession}}.
  \describe{
    \item{\code{\link{browserView}(object, range =
        range(object), track = trackNames(object), ...)}}{
      Constructs a \code{\linkS4class{BrowserView}} of \code{range}
      for this session.
    }
    \item{\code{\link{browserViews}(object, ...)}}{
      Gets the \code{\linkS4class{BrowserView}} instances belonging to
      this session.}
    \item{\code{\link{activeView}(object, ...)}}{
      Returns the \code{\linkS4class{BrowserView}} that is currently
      active in the session.}
    \item{\code{\link{range}(x, ...)}}{
      Gets the \code{\link[IRanges:RangesList-class]{RangesList}}
      representing the range of the genome currently displayed by the
      browser (i.e. the range shown by the active view) or a default
      value (possibly \code{NULL}) if no views exist.
    }
    \item{\code{\link[BSgenome]{getSeq}(object, range =
        range(object), ...)}}{
      gets a genomic sequence of \code{range} from this session.}
    \item{\code{\link[=sequence<-]{sequence}(object, ...) <- value}}{
      Loads a sequence into the session.}
    \item{\code{\link[=track<-]{track}(object, name =
        deparse(substitute(track)), view = TRUE, ...) <- value}}{
      Loads one or more \code{track}s into the session and optionally open a
      \code{view} of the track.
    }
    \item{\code{x[[i]] <- value}}{Loads the track \code{value} into
      session \code{x}, under the name \code{i}. Shortcut to above.}
    \item{\code{x$name <- value}}{Loads the track \code{value} into
      session \code{x}, under the name \code{name}. Shortcut to above.}
    \item{\code{\link{track}(object, ...)}}{
      Gets a track from a session as a
      \code{\link[IRanges:RangedData-class]{RangedData}}.
    }
    \item{\code{x[[i]]}}{Gets the track named \code{i} from session
      \code{x}. A shortcut to \code{track}.
    }
    \item{\code{x$name}}{Gets the track named \code{name} from session
      \code{x}. A shortcut to \code{track}.}
    \item{\code{\link{trackNames}(object, ...)}}{ Gets the names of the
      tracks stored in this session. }
    \item{\code{\link{genome}(x)}, \code{genome(x) <- value}}{Gets or
      sets the genome identifier (e.g. \dQuote{hg18}) for the session.}
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
