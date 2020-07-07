\name{FastaFile-class}
\docType{class}

%% Classes:
\alias{class:FastaFile}
\alias{FastaFile-class}

%% Constructor:
\alias{FastaFile}

%% Import:
\alias{import,FastaFile,ANY,ANY-method}

%% Export:
\alias{export,ANY,FastaFile,ANY-method}
\alias{export,XStringSet,FastaFile,ANY-method}

\title{FastaFile objects}

\description{
  These functions support the import and export of the Fasta
  sequence format, using the Biostrings package.
}

\usage{
\S4method{import}{FastaFile,ANY,ANY}(con, format, text,
           type = c("DNA", "RNA", "AA", "B"), ...)

\S4method{export}{ANY,FastaFile,ANY}(object, con, format, ...)
\S4method{export}{XStringSet,FastaFile,ANY}(object, con, format, ...)
}

\arguments{
  \item{con}{A path or \code{FastaFile} object.  URLs and connections
    are not supported. If \code{con} is not a \code{FastaFile}, either
    the file extension or the \code{format} argument needs to be
    \dQuote{fasta}. Compressed files (\dQuote{gz}, \dQuote{bz2} and
    \dQuote{xz}) are handled transparently.
  }
  \item{object}{The object to export, should be an \code{XStringSet} or
    something coercible to a \code{DNAStringSet}, like a character vector.
  }
  \item{format}{If not missing, should be \dQuote{fasta}.
  }
  \item{text}{If \code{con} is missing, a character vector to use as the
    input
  }
  \item{type}{Type of biological sequence.}
  \item{...}{Arguments to pass down
    to \code{\link[Biostrings]{writeXStringSet}} (export)
    or the \code{\link[Biostrings]{readDNAStringSet}} family of
    functions (import).
  }
}

\seealso{
  These functions are implemented by the Biostrings
  \code{\link[Biostrings]{writeXStringSet}} (export) and the
  \code{\link[Biostrings]{readDNAStringSet}} family of functions
  (import).

  See \link[BSgenome]{export-methods} in the \pkg{BSgenome} package for
  exporting a \link[BSgenome]{BSgenome} object as a FASTA file.
}

\section{FastaFile objects}{
  The \code{FastaFile} class extends
  \code{\link[BiocIO:BiocFile-class]{BiocFile}} and is
  a formal represention of a resource in the Fasta format.
  To cast a path, URL or connection to a \code{FastaFile}, pass it to
  the \code{FastaFile} constructor.
}

\author{Michael Lawrence}

\keyword{methods}
\keyword{classes}
