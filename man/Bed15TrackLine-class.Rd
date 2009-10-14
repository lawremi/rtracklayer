\name{Bed15TrackLine-class}
\docType{class}
\alias{Bed15TrackLine-class}
\alias{coerce,character,Bed15TrackLine-method}
\alias{coerce,Bed15TrackLine,character-method}

\title{Class "Bed15TrackLine"}
\description{ A \acronym{UCSC} track line for graphical tracks. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Bed15TrackLine",
    ...)} or parsed from a character vector track line with
  \code{as(text, "Bed15TrackLine")}.
}
\section{Slots}{
  \describe{
    \item{\code{expStep}:}{A \code{"numeric"} scalar indicating the step
      size for the heatmap color gradient.}
    \item{\code{expScale}:}{A positive \code{"numeric"} scalar indicating the
      range of the data to be \code{[-expScale, expScale]} for
      determining the heatmap color gradient.}
    \item{\code{expNames}:}{A \code{"character"} vector naming the
      the experimental samples.}
    \item{\code{name}:}{Object of class \code{"character"} specifying
      the name of the track. }
    \item{\code{description}:}{Object of class \code{"character"}
      describing the track. }
    \item{\code{visibility}:}{Object of class \code{"character"}
      indicating the default visible mode of the track, see
      \code{\linkS4class{UCSCTrackModes}}.}
    \item{\code{color}:}{Object of class \code{"integer"} representing the
      track color (as from \code{\link{col2rgb}}).
    }
    \item{\code{priority}:}{Object of class \code{"numeric"} specifying
      the rank of this track. }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{TrackLine}"}, directly.
}
\section{Methods}{
  \describe{
    \item{as(object, "character")}{Export line to its string
      representation.}
  }
}
\references{ Official documentation:
  \url{http://genomewiki.ucsc.edu/index.php/Microarray_track}. }
\author{ Michael Lawrence }
\seealso{
  \code{\link{export.bed15}} for exporting bed15 tracks.
}
\keyword{classes}
