\name{ucscSession-class}
\docType{class}
\alias{ucscSession-class}
\alias{initialize,ucscSession-method}

\title{Class "ucscSession"}
\description{An implementation of \code{\linkS4class{browserSession}}
  for the \acronym{UCSC} genome browser.}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{\link{browserSession}("ucsc", url =
    "http://genome.ucsc.edu/cgi-bin")}.
}
\section{Slots}{
  \describe{
    \item{\code{url}:}{Object of class \code{"character"} holding the
      base \acronym{URL} of the \acronym{UCSC} browser. }
    \item{\code{hguid}:}{Object of class \code{"numeric"} holding the
      user identification code.}
    \item{\code{views}:}{Object of class \code{"environment"} containing
    a list stored under the name "instances". The list holds the instances of
    \code{\linkS4class{browserView}} for this session. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{browserSession}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{browserView}(object, segment =
        genomeSegment(object), track = tracks(object), ...)}}{ Creates a
      \code{\linkS4class{browserView}} of \code{segment} with visible
      tracks specified by \code{track}. \code{track} may be an instance
      of \code{\linkS4class{ucscTrackModes}}. Arguments in \code{...} should
      override slots in \code{segment} or else match parameters
      to a \code{\link{ucscTrackModes}} method for creating a
      \code{ucscTrackModes} instance that will override modes indicated
      by the \code{track} parameter. }
    \item{\code{\link{browserViews}}(object)}{ Gets the
      \code{\linkS4class{browserView}} instances for this session. }
    \item{\code{\link{genomeSegment}}(object)}{ Gets the
      \code{\linkS4class{genomeSegment}} last displayed in this session. }
    \item{\code{\link{genomeSequence}(object, segment, track = "Assembly")}}{
      Gets the sequence in \code{segment} and \code{track}.
    }
    \item{\code{\link{layTrack}(object, track, name = names(track),
        view = TRUE, format = "gff", ...)}}{
      Loads a track, stored under \code{name} and formatted as
      \code{format}. The arguments in \code{...} are passed on to
      \code{\link{export.ucsc}}, so they could be slots in a
      \code{\linkS4class{ucscTrackLine}} subclass or parameters to pass
      on to the export function for \code{format}.
    }
    \item{\code{\link{trackSet}(object, segment = genomeSegment(), name =
        deparse(substitute(object)), table = NULL)}}{
      Retrieves a \code{\linkS4class{trackSet}} with features in
      \code{segment} from track named \code{name}. Some built-in tracks
      have multiple series, each stored in a separate database
      table. A specific table may be retrieved by passing its name in
      the \code{table} parameter.
    }
    \item{\code{\link{tracks}}(object)}{ Gets the
      names of the tracks stored in the session. }
    \item{\code{\link{ucscTrackModes}}(object)}{ Gets the
      default view modes for the tracks in the session. }
    \item{\code{\link{ucscTable}}(object)}{ Gets the database \code{table}
      in \code{segment} from \code{track} as a \code{data.frame}. }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{browserSession}} for creating instances of this class.
}
\keyword{classes}
