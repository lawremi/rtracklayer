\name{RTLFile-class}
\docType{class}

%% Classes:
\alias{class:RTLFile}
\alias{RTLFile-class}
\alias{class:CompressedFile}
\alias{CompressedFile-class}
\alias{RTLFileList-class}

%% Accessors
\alias{path}
\alias{path,RTLFile-method}
\alias{resource}
\alias{fileFormat}
\alias{fileFormat,RTLFile-method}
\alias{fileFormat,CompressedFile-method}
\alias{fileFormat,character-method}

%% Show:
\alias{show,RTLFile-method}

%% Related functions:
\alias{FileForFormat}

%% Coerce:
\alias{as.character,RTLFile-method}

\title{RTLFile objects}

\description{
  A \code{RTLFile} object is the base class for classes representing
  files accessible with rtracklayer. It wraps a resource (either a path,
  URL or connection). We can represent a list of \code{RTLFile} objects
  with a \code{RTLFileList}.
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
    \item{}{
      \code{fileFormat(x)}: Gets a string identifying the file
      format. Can also be called directly on a character file path, in
      which case it uses a heuristic based on the file extension.
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{
      \code{as.character(x)}:
      Returns the path of the file as a character vector.
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
