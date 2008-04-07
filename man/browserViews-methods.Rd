\name{browserViews-methods}
\docType{methods}
\alias{browserViews}
\alias{browserViews-methods}
\alias{browserViews,ucscSession-method}
\alias{browserViews,argoSession-method}
\title{ Getting the browser views }
\description{
 Methods for getting browser views.
}
\section{Methods}{
\describe{
  The following methods are defined by \pkg{rtracklayer}. 
  \item{object = "ucscSession"}{ Gets the instances of
    \code{\linkS4class{browserView}} in the session. }

  \item{object = "argoSession"}{ Gets the instances of
    \code{\linkS4class{browserView}} in the session. }
}}
\seealso{\code{\link{browserView}} for creating a browser view.}
\examples{
\dontrun{
session <- browseGenome()
browserViews(session)
}
}
\keyword{methods}
