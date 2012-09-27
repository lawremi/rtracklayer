\name{TrackDb-class}
\docType{class}

%% Classes:
\alias{class:TrackDb}
\alias{TrackDb-class}

%% Accessors:
\alias{$,TrackDb-method}
\alias{$<-,TrackDb-method}
\alias{[[,TrackDb,ANY,ANY-method}
\alias{[[<-,TrackDb,ANY,ANY-method}

%% Data Access:
\alias{track<-,TrackDb,ANY-method}
\alias{track<-,TrackDb,RangedData-method}
\alias{track<-,TrackDb,RangedDataList-method}

\title{Track Databases}

\description{
  The \code{TrackDb} class is an abstraction around a database of
  tracks. Implementations include \code{\linkS4class{BrowserSession}}
  derivatives and \code{\linkS4class{QuickloadGenome}}. Here, a track is
  defined as an interval dataset.
}

\section{Accessor Methods}{
  Every implementation should support these methods:
  
  \describe{
    \item{}{
      \code{length(x)}: number of tracks
    }
    \item{}{
      \code{names(x)}, \code{trackNames(x)}: names of the tracks
    }
    \item{}{
      \code{mcols(x)}: merged metadata on the tracks
    }
    \item{}{
      \code{track(x, name)}, \code{x$name}, \code{x[[name]]}:
      get the track called \code{name}
    }
    \item{}{
      \code{track(x, name) <- value}, \code{x$name <- value},
      \code{x[[name]] <- value}: store the track \code{value} under
      \code{name}. Different implementations will support different
      types for \code{value}. Generally, an interval data structure like
      \code{GenomicRanges}.
    }
  }
}

\author{Michael Lawrence}

\keyword{methods}
\keyword{classes}
