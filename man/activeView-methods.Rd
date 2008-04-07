\name{activeView-methods}
\docType{methods}
\alias{activeView}
\alias{activeView<-}
\alias{activeView-methods}
\alias{activeView<--methods}
\alias{activeView,browserSession-method}
\alias{activeView,argoSession-method}
\alias{activeView<-,argoSession-method}
\alias{activeView,argoView-method}
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
    \item{object = "ucscView" or "argoView"}{
      \code{activeView(object)}: Gets a logical indicating whether this
      is the active view.
    }
}}
\examples{
\dontrun{
  session <- browseGenome(browser = "argo")
  activeView(session)
  view <- browserView(session, genomeSegment("hg18", "chr22", 20000, 30000))
  activeView(session) <- view
}}
\keyword{methods}
