\name{argoView-class}
\docType{class}
\alias{argoView-class}
\alias{close,argoView-method}

\title{Class "argoView"}
\description{An object representing a view of a genome in the Argo browser.}
\section{Objects from the Class}{
  Objects may be created through calls of the form
  \code{\link{browserView}(session, ...)} where \code{session} is an
  instance of \code{\linkS4class{argoSession}}.
}
\section{Slots}{
  \describe{
    \item{\code{jref}:}{Object of class \code{"jobjRef"} holding a
      reference to the low-level Argo object. }
    \item{\code{session}:}{Object of class \code{"browserSession"} to
      which this view belongs. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{browserView}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\link{activeView}(object)}}{ Returns a logical
      indicating whether this view is the active view. }
    \item{\code{\link{close}(con)}}{ Closes this view. }
    \item{\code{\link{genomeSegment}(object)}}{ Obtains the
      \code{\linkS4class{genomeSegment}} displayed by this view. }
    \item{\code{\link{genomeSegment}(object) <- value}}{ Sets the
      \code{\linkS4class{genomeSegment}} displayed by the view. }
    \item{\code{\link{tracks}(object)}}{ Gets the names of the visible tracks. }
    \item{\code{\link{tracks<-}(object)}}{ Sets the visible tracks by
      name. Note that in Argo track visibility is global. }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\linkS4class{browserView}} for creating instances from a
  \code{\linkS4class{argoSession}}.
}
\keyword{classes}
