\name{RTLFile-class}
\docType{class}

%% Classes:
\alias{class:RTLFile}
\alias{RTLFile-class}
\alias{class:CompressedFile}
\alias{CompressedFile-class}

%% Accessors
\alias{path}
\alias{path,RTLFile-method}
\alias{resource}

%% Show:
\alias{show,RTLFile-method}

%% Related functions:
\alias{FileForFormat}

\title{RTLFile objects}

\description{
  A \code{RTLFile} object is the base class for classes representing
  files accessible with rtracklayer. It wraps a resource
  (either a path, URL or connection).
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{RTLFile}
  object.
  
  \describe{
    \item{}{
      \code{path(x)}:
      Gets the path, as a \code{character} vector, to the resource
      represented by the \code{RTLFile} object, if possible.
    }
    \item{}{
      \code{resource(x)}:
      Gets the low-level resource, either a character vector (a path or
      URL) or a connection.
    }
  }
}

\section{Related functions}{
  \describe{
    \item{}{
      \code{FileForFormat(path, format = file_ext(path))}:
      Determines the file type of \code{path} and returns
      a high-level file object such as BamFile, BEDFile,
      BigWigFile etc.. 
    }
  }
}

\author{Michael Lawrence}

\seealso{
  Implementing classes include: \code{\linkS4class{BigWigFile}},
  \code{\linkS4class{TwoBitFile}}, \code{\linkS4class{BEDFile}},
  \code{\linkS4class{GFFFile}}, and \code{\linkS4class{WIGFile}}.
}

\keyword{methods}
\keyword{classes}
