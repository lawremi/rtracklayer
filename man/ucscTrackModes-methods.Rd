\name{ucscTrackModes-methods}
\docType{methods}
\alias{ucscTrackModes}
\alias{ucscTrackModes<-}
\alias{ucscTrackModes-methods}
\alias{ucscTrackModes,character-method}
\alias{ucscTrackModes,missing-method}
\alias{ucscTrackModes,UCSCView-method}
\alias{ucscTrackModes<--methods}
\alias{ucscTrackModes<-,UCSCView,UCSCTrackModes-method}
\alias{ucscTrackModes<-,UCSCView,character-method}
\alias{ucscTrackModes,ucscTracks-method}
\alias{ucscTrackModes,UCSCSession-method}

\title{ Accessing UCSC track modes }
\description{
  Generics for getting and setting UCSC track visibility modes ("hide",
  "dense", "full", "pack", "squish").
}
\section{Methods}{
\describe{
  The following methods are defined by \pkg{rtracklayer} for \bold{getting} the
  track modes through the generic \code{ucscTrackModes(object, ...)}.
  
  \item{object = "character"}{
    \code{function(object, hide = character(),
                   dense = character(), pack = character(),
                   squish = character(), full = character())}
    Creates an instance of \code{\linkS4class{UCSCTrackModes}} from
    \code{object}, a character vector of mode names, with the
    corresponding track ids given in the \code{names} attribute. Note
    that \code{object} can be a \code{UCSCTrackModes} instance, as
    \code{UCSCTrackModes} extends \code{character}.
    The other parameters are character vectors identifying the tracks for
    each mode and overriding the modes specified by
    \code{object}.
  }

  \item{object = "missing"}{
    The same interface as above, except \code{object} defaults to an empty
    character vector.
  }

  \item{object = "UCSCView"}{
    Gets modes for tracks in the view.
  }
  
  \item{object = "UCSCSession"}{
    Gets default modes for the tracks in the session. These are the
    modes that will be used as the default for a newly created view.
  }

  The following methods are defined by \pkg{rtracklayer} for \bold{setting} the
  track modes through the generic \code{ucscTrackModes(object) <- value}.

  \item{object = "UCSCView", value = "UCSCTrackModes"}{
    Sets the modes for the tracks in the view.
  }

  \item{object = "UCSCView", value = "character"}{
    Sets the modes from a character vector of mode names, with the
    corresponding track names given in the \code{names} attribute.
  }
}}

\seealso{\code{\link{trackNames}} and \code{\link{trackNames<-}} for just
  getting or setting which tracks are visible (not of mode "hide").}

\examples{
# Tracks "foo" and "bar" are fully shown, "baz" is hidden
modes <- ucscTrackModes(full = c("foo", "bar"), hide = "baz")
# Update the modes to hide track "bar"
modes2 <- ucscTrackModes(modes, hide = "bar")
}
\keyword{methods}
