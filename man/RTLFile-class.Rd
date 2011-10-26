\name{RTLFile-class}
\docType{class}

%% Classes:
\alias{class:RTLFile}
\alias{RTLFile-class}

%% Accessors
\alias{path}
\alias{path,RTLFile-method}

%% Show:
\alias{show,RTLFile-method}

\title{RTLFile objects}

\description{
  A \code{RTLFile} object is the base class for classes representing
  files accessible with rtracklayer. It currently stores a path and
  provides a few utilities.
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{RTLFile}
  object.
  
  \describe{
    \item{}{
      \code{path(x)}:
      Gets the path, as a \code{character} vector, to the file
      represented by the \code{RTLFile} object.
    }
  }
}


\author{Michael Lawrence}

\seealso{
  Implementing classes like \code{\linkS4class{BigWigFile}} and
  \code{\linkS4class{TwoBitFile}}
}

\keyword{methods}
\keyword{classes}
