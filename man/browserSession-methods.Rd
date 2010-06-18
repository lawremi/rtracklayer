\name{browserSession-methods}
\docType{methods}
\alias{browserSession}
\alias{browserSession-methods}
\alias{browserSession,character-method}
\alias{browserSession,BrowserView-method}
\alias{browserSession,missing-method}
\title{ Get a genome browser session }
\description{
  Methods for getting browser sessions.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
  
  \describe{
    \item{object = "character"}{ \code{browserSession(object, ...)}: Creates a
      \code{\linkS4class{BrowserSession}} from a genome browser
      identifier. The identifier corresponds to the prefix of the
      session class name (e.g. "UCSC" in "UCSCSession"). The arguments
      in \dots are passed to the initialization function of the class.}
    
    \item{object = "browserView"}{ Gets the
      \code{\linkS4class{BrowserSession}} for the view. }

    \item{object = "missing"}{ Calls \code{browserSession("ucsc", ...)}. }
  }
}
\keyword{methods}
