\name{BrowserView-class}
\docType{class}
\alias{BrowserView-class}
\alias{show,BrowserView-method}

\title{Class "BrowserView"}
\description{ An object representing a genome browser view of a particular
  \code{\linkS4class{RangesList}} of a genome. }
\section{Objects from the Class}{
  A virtual Class: No objects may be created from it directly. See
  \code{\link{browserView}} for obtaining
  an instance of an implementation for a particular genome browser.
}
\section{Slots}{
  \describe{
    \item{\code{session}:}{Object of class \code{"BrowserSession"} the
      browser session to which this view belongs. }
  }
}
\section{Methods}{
  This specifies the API implemented by each browser backend. Note that
  a backend is not guaranteed to support all operations.
  \describe{
    \item{\code{\link{browserSession}(object)}}{ Obtains the
      \code{\linkS4class{BrowserSession}} to which this view belongs. }
    \item{\code{\link{close}(object)}}{ Close this view. }
    \item{\code{\link{range}(object)}}{ Obtains the
      \code{\linkS4class{RangesList}} displayed by this view.}
    \item{\code{\link{trackNames}(object)}}{ Gets the names of the visible
      tracks in the view. }
    \item{\code{\link{trackNames}(object) <- value}}{ Sets the visible
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
