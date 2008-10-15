\name{browserSession-methods}
\docType{methods}
\alias{browserSession}
\alias{browserSession-methods}
\alias{browserSession,character-method}
\alias{browserSession,browserView-method}
\alias{browserSession,missing-method}
\title{ Get a genome browser session }
\description{
  Methods for getting browser sessions.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}. 
  
  \describe{
    \item{object = "character"}{ \code{browserSession(object, ...)}: Creates a
      \code{\linkS4class{browserSession}} from a genome browser
      identifier. The identifier corresponds to the prefix of the
      session class name (e.g. "ucsc" in "ucscSession"). The arguments
      in \dots are passed to the initialization function of the class.}
    
    \item{object = "browserView"}{ Gets the
      \code{\linkS4class{browserSession}} for the view. }

    \item{object = "missing"}{ Calls \code{browserSession("ucsc", ...)}. }
  }
}
\keyword{methods}
