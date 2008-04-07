\name{ucscTrackModes-class}
\docType{class}
\alias{ucscTrackModes-class}

\title{Class "ucscTrackModes"}
\description{A vector of view modes ("hide", "dense", "full",
  "pack", "squish") for each track in a \acronym{UCSC} view.}
\section{Objects from the Class}{
  Objects may be created by calls of the form
  \code{\link{ucscTrackModes}(object =
    character(), hide = character(), dense = character(), pack = character(),
    squish = character(), full = character())},
  where \code{object} should be a character vector of mode names (with
  its \code{names} attribute specifying the corresponding track names).
  The other parameters should contain track names that override the
  modes in \code{object}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} holding the
      modes ("hide", "dense", "full", "pack", "squish"), with its
      \code{names} attribute holding corresponding track names.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
Class \code{"\link[Biobase:characterORMIAME-class]{characterORMIAME}"}, by class "character", distance 2.
}
\section{Methods}{
  \describe{
    \item{\code{\link{tracks}(object)}}{ Gets the names of the visible
      tracks (those that do not have mode "hide"). }
    \item{\code{\link{tracks}(object) <- value}}{ Sets the names of the
      visible tracks. Any tracks named in \code{value} are set to "full"
      if the are currently set to "hide" in this object. Any tracks not
      in \code{value} are set to "hide". All other modes are preserved.}
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\linkS4class{ucscView}} on which track view modes may be set.
}
\keyword{classes}
