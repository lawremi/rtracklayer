\name{Quickload-class}
\docType{class}

%% Classes:
\alias{class:Quickload}
\alias{Quickload-class}

%% Constructor:
\alias{Quickload}

%% Accessors:
\alias{$,Quickload-method}
\alias{[[,Quickload,ANY,ANY-method}
\alias{genome,Quickload-method}
\alias{length,Quickload-method}
\alias{uri}

%% Coercion:
\alias{coerce,character,Quickload-method}

%% Show:
\alias{show,Quickload-method}

\title{Quickload Access}

\description{
  The \code{Quickload} class represents a Quickload data source,
  essentially directory layout separating tracks and sequences by
  genome, along with a few metadata files. This interface abstracts
  those details and provides access to a Quickload at any URL supported
  by R (HTTP, FTP, and local files). This is an easy way to make data
  accessible to the Integrated Genome Browser (IGB).
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{Quickload(uri = "quickload", create = FALSE)}: Constructs a
      new \code{Quickload} object, representing a repository at
      \code{uri}. If \code{create} is \code{TRUE}, and \code{uri} is
      writeable (i.e., local), the repository is created if it does not
      already exist. If it does exist, then a message is emitted to
      indicate that the repository was not recreated.
    }
  }
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{Quickload}
  object.
  
  \describe{
    \item{}{\code{x$genome}, \code{x[["genome"]]}: Get
      the \code{\linkS4class{QuickloadGenome}} object for the genome named
      \code{genome}. This is where all the data is stored.
    }
    \item{}{\code{length(x)}: number of genomes in the repository
    }
    \item{}{\code{uri(x)}:
      Get the URI pointing to the Quickload repository.
    }
    \item{}{\code{genome(x)}:
      Get the identifiers of the genomes present in the repository.
    }
  }
}

\author{Michael Lawrence}

\examples{
ql <- Quickload(system.file("tests", "quickload", package = "rtracklayer"))
uri(ql)
genome(ql)
ql$T_species_Oct_2011
}

\keyword{methods}
\keyword{classes}
