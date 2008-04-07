\name{browserView-class}
\docType{class}
\alias{browserView-class}
\alias{show,browserView-method}

\title{Class "browserView"}
\description{ An object representing a genome browser view of a particular
  \code{\linkS4class{genomeSegment}} of a genome. }
\section{Objects from the Class}{
  A virtual Class: No objects may be created from it directly. See
  \code{\link{browserView}} for obtaining
  an instance of an implementation for a particular genome browser.
}
\section{Slots}{
  \describe{
    \item{\code{session}:}{Object of class \code{"browserSession"} the
      browser session to which this view belongs. }
  }
}
\section{Methods}{
  This specifies the API implemented by each browser backend. Note that
  a backend is not guaranteed to support all operations.
  \describe{
    \item{\code{\link{browserSession}(object)}}{ Obtains the
      \code{\linkS4class{browserSession}} to which this view belongs. }
    \item{\code{\link{close}(object)}}{ Close this view. }
    \item{\code{\link{genomeSegment}(object)}}{ Obtains the
      \code{\linkS4class{genomeSegment}} displayed by this view.}
    \item{\code{\link{tracks}(object)}}{ Gets the names of the visible
      tracks in the view. }
    \item{\code{\link{tracks}(object) <- value}}{ Sets the visible
      tracks by their names. }
    \item{\code{\link{show}(object)}}{ Outputs a textual description of
      this view. }  
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{browserView}} for obtaining instances of this class.
}
\keyword{classes}
