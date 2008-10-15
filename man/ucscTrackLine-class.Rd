\name{ucscTrackLine-class}
\docType{class}
\alias{ucscTrackLine-class}
\alias{coerce,character,ucscTrackLine-method}
\alias{coerce,ucscTrackLine,character-method}
\alias{show,ucscTrackLine-method}

\title{Class "ucscTrackLine"}
\description{An object representing a "track line" in the \acronym{UCSC}
  format. There are two concrete types of track lines:
  \code{\linkS4class{basicTrackLine}} (used for most types of tracks) and
  \code{\linkS4class{wigTrackLine}} (used for Wiggle tracks). This class
  only declares the common elements between the two.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ucscTrackLine",
    ...)} or parsed from a character vector track line with
  \code{as(text, "ucscTrackLine")}.
  But note that \acronym{UCSC} only understands one of the
  subclasses mentioned above.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} specifying
      the name of the track. }
    \item{\code{description}:}{Object of class \code{"character"}
      describing the track. }
    \item{\code{visibility}:}{Object of class \code{"character"}
      indicating the default visible mode of the track, see
      \code{\linkS4class{ucscTrackModes}}.}
    \item{\code{color}:}{Object of class \code{"integer"} representing the
      track color (as from \code{\link{col2rgb}}).
    }
    \item{\code{priority}:}{Object of class \code{"numeric"} specifying
      the rank of this track. }
  }
}
\section{Methods}{
  \item{as(object, "character")}{Export line to its string representation.}
}
\references{
  \url{http://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK}
  for the official documentation. }
\author{ Michael Lawrence }

\seealso{
  \code{\linkS4class{basicTrackLine}} (used for most types of tracks) and
  \code{\linkS4class{wigTrackLine}} (used for Wiggle tracks).
}
\keyword{classes}
