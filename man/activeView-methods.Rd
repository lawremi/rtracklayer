\name{activeView-methods}
\docType{methods}
\alias{activeView}
\alias{activeView<-}
\alias{activeView-methods}
\alias{activeView<--methods}
\alias{activeView,browserSession-method}
\alias{activeView,ucscView-method}
\title{Accessing the active view}
\description{
 Get the active view.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}.
  \describe{
    \item{object = "browserSession"}{
      \code{activeView(object)}: Gets the active
      \code{\linkS4class{browserView}} from a browser session.

      \code{activeView(object) <- value}: Sets the active
      \code{\linkS4class{browserView}} in a browser session.
    }
}}
\keyword{methods}
