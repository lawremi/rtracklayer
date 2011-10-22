\name{UCSCSession-class}
\docType{class}
\alias{UCSCSession-class}
\alias{initialize,UCSCSession-method}
\alias{track,UCSCSession-method}
\alias{range<-}
\alias{range<-,UCSCSession-method}
\alias{genome,UCSCSession-method}
\alias{genome<-,UCSCSession-method}
\alias{range,UCSCSession-method}
\alias{seqinfo,UCSCSession-method}
\alias{range,ucscCart-method}
\alias{genome,ucscCart-method}

\title{Class "UCSCSession"}
\description{An implementation of \code{\linkS4class{BrowserSession}}
  for the \acronym{UCSC} genome browser.}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{\link{browserSession}("ucsc", url =
    "http://genome.ucsc.edu/cgi-bin", ...)}. The arguments in \code{...}
  correspond to libcurl options, see
  \code{\link[RCurl]{getCurlHandle}}. Setting these options may be
  useful e.g. for getting past a proxy.
}
\section{Slots}{
  \describe{
    \item{\code{url}:}{Object of class \code{"character"} holding the
      base \acronym{URL} of the \acronym{UCSC} browser. }
    \item{\code{hguid}:}{Object of class \code{"numeric"} holding the
      user identification code.}
    \item{\code{views}:}{Object of class \code{"environment"} containing
    a list stored under the name "instances". The list holds the instances of
    \code{\linkS4class{BrowserView}} for this session. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{BrowserSession}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{browserView}(object, range =
        range(object), track = trackNames(object), ...)}}{ Creates a
      \code{\linkS4class{BrowserView}} of \code{range} with visible
      tracks specified by \code{track}. \code{track} may be an instance
      of \code{\linkS4class{UCSCTrackModes}}. Arguments in \code{...} should
      match parameters
      to a \code{\link{ucscTrackModes}} method for creating a
      \code{UCSCTrackModes} instance that will be merged with and
        override modes indicated by the \code{track} parameter. }
    \item{\code{\link{browserViews}}(object)}{ Gets the
      \code{\linkS4class{BrowserView}} instances for this session. }
    \item{\code{\link{range}}(x, asRangedData = TRUE)}{ Gets the
      \code{\link[GenomicRanges:GRanges-class]{GRanges}} last
      displayed in this session. Set \code{asRangedData}
      to \code{FALSE} to obtain a \code{GRanges} object.}
    \item{\code{\link[rtracklayer:RangedData-methods]{genome}(x)}}{Gets
      the genome identifier of the session,
      i.e. \code{genome(range(x))}.}
    \item{\code{seqinfo}}{ Gets
      the \code{\link[GenomicRanges:Seqinfo]{Seqinfo}} object with the
      lengths of the chromosomes in the currenet genome. No
      circularity information is available.
    }
    \item{\code{range(x) <- value}}{ Sets \code{value}, usually a
      \code{GRanges} object or \code{RangesList}, as the range
      of session \code{x}. Note that
      this setting only lasts until a view is created or
      manipulated. This mechanism is useful, for example, when treating the UCSC
      browser as a database, rather than a genome viewer.}
    \item{\code{\link[rtracklayer:RangedData-methods]{genome}(x) <- value}}{ 
      Sets the genome identifier on the range of session \code{x}.}
    \item{\code{\link[Biostrings]{getSeq}(object, range, track = "Assembly")}}{
      Gets the sequence in \code{range} and \code{track}.
    }
    \item{\code{\link[=track<-]{track}(object, name = names(track),
        format = "auto", ...) <- value}}{
      Loads a track, stored under \code{name} and formatted as
      \code{format}. The "auto" format resolves to "bed" for qualitative
        data. For quantitative data, i.e., data with a numeric \code{score}
        column, "wig" or "bedGraph" is chosen, depending on how well the
        data compresses into wig. The arguments in \code{...} are passed on to
      \code{\link{export.ucsc}}, so they could be slots in a
      \code{\linkS4class{TrackLine}} subclass (and thus specify
        visual attributes like color) or parameters to pass
      on to the export function for \code{format}.
    }
    \item{\code{\link{track}(object, name, range = range(object),
        table = NULL)}}{
      Retrieves a \code{\link[IRanges:RangedData-class]{RangedData}}
      with features in 
      \code{range} from track named \code{name}. Some built-in tracks
      have multiple series, each stored in a separate database
      table. A specific table may be retrieved by passing its name in
      the \code{table} parameter. See \code{\link{tableNames}} for a way
      to list the available tables.
    }
    \item{\code{\link{trackNames}}(object)}{ Gets the
      names of the tracks stored in the session. }
    \item{\code{\link{ucscTrackModes}}(object)}{ Gets the
      default view modes for the tracks in the session. }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{browserSession}} for creating instances of this class.
}
\keyword{classes}
