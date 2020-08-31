\name{TrackHub-class}
\docType{class}

%% Classes:
\alias{class:TrackHub}
\alias{TrackHub-class}

%% Constructor:
\alias{TrackHub}

%% Accessors:
\alias{$,TrackHub-method}
\alias{[[,TrackHub,ANY,ANY-method}
\alias{genome,TrackHub-method}
\alias{length,TrackHub-method}
\alias{uri}
\alias{writeTrackHub}

%% Data Access:
\alias{hub,TrackHub-method}
\alias{hub<-,TrackHub-method}
\alias{shortLabel,TrackHub-method}
\alias{shortLabel<-,TrackHub-method}
\alias{longLabel,TrackHub-method}
\alias{longLabel<-,TrackHub-method}
\alias{genomeFile,TrackHub-method}
\alias{genomeFile<-,TrackHub-method}
\alias{email,TrackHub-method}
\alias{email<-,TrackHub-method}
\alias{descriptionUrl,TrackHub-method}
\alias{descriptionUrl<-,TrackHub-method}
\alias{genomeField,TrackHub-method}
\alias{genomeField<-,TrackHub-method}
\alias{genomeInfo,TrackHub-method}
\alias{genomeInfo<-,TrackHub-method}

%% Coercion:
\alias{coerce,character,TrackHub-method}

%% Show:
\alias{show,TrackHub-method}

\title{TrackHub Access}

\description{
  The \code{TrackHub} class represents a TrackHub data source,
  essentially directory layout separating tracks and sequences by
  genome, along with a few metadata files. This interface abstracts
  those details and provides access to a TrackHub at any URL supported
  by R (HTTP, FTP, and local files). This is an easy way to make data
  accessible to the UCSC Genome Browser.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{TrackHub(uri, create = FALSE)}: Constructs a
      new \code{TrackHub} object, representing a repository at
      \code{uri}. If \code{create} is \code{TRUE}, and \code{uri} is
      writeable (i.e., local), the repository is created if it does not
      already exist. If it does exist, then a message is emitted to
      indicate that the repository was not recreated.
    }
  }
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{TrackHub}
  object.

  \describe{
    \item{}{\code{x$genome}, \code{x[["genome"]]}: Get
      the \code{\linkS4class{TrackHubGenome}} object for the genome named
      \code{genome}.
    }
    \item{}{\code{length(x)}: number of genomes in the repository.
    }
    \item{}{\code{uri(x)}:
      Get the URI pointing to the TrackHub repository.
    }
    \item{}{\code{genome(x)}:
      Get the identifiers of the genomes present in the repository.
    }
    \item{}{\code{writeTrackHub(x)}:
      Write hub content and genomes from memory representation to the hub file and genomes file.
      It also create resources if they are missing like genomes file and genome directory for
      newly add genome.
    }
  }
}

\section{Data Access}{
  Note that all storing methods(like \code{hub()<-}) are only supported for local repositories, i.e., those with a file:// URI scheme.
  \describe{
    \item{}{
      \code{hub(x)}: get the value of hub.
    }
    \item{}{
      \code{hub(x) <- value}: store the \code{value} of hub for \code{x}.
    }
    \item{}{
      \code{shortLabel(x)}: get the value of hub.
    }
    \item{}{
      \code{shortLabel(x) <- value}: store the \code{value} of shortLabel for \code{x}.
    }
    \item{}{
      \code{longLabel(x)}: get the value of hub.
    }
    \item{}{
      \code{longLabel(x) <- value}: store the \code{value} of longLabel for \code{x}.
    }
    \item{}{
      \code{genomeFile(x)}: get the value of hub.
    }
    \item{}{
      \code{genomeFile(x) <- value}: store the \code{value} of genomesFile for \code{x}.
    }
    \item{}{
      \code{email(x)}: get the value of hub.
    }
    \item{}{
      \code{email(x) <- value}: store the \code{value} of email for \code{x}.
    }
    \item{}{
      \code{descriptionUrl(x)}: get the value of hub.
    }
    \item{}{
      \code{descriptionUrl(x) <- value}: store the \code{value} of descriptionUrl for \code{x}.
    }
    \item{}{
      \code{genomeField(x, name, field)}: Get the \code{value} of \code{field} for \code{name} genome.
    }
    \item{}{
      \code{genomeField(x, name, field) <- value}: Set or Update the \code{field} and \code{value} for \code{name} genome.
    }
    \item{}{
      \code{genomeInfo(x, name)}: Get the \code{Genome} object for \code{name} genome.
    }
    \item{}{
      \code{genomeInfo(x) <- value}: Add \code{value} (Genome object) to existing genomes list.
      \code{Genome} takes named arguemnts of all UCSC supported fields for genome
      file(like \code{genome, trackDb, twoBitPath}, etc).
    }
  }
}

\author{Michael Lawrence}

\examples{
th <- TrackHub(system.file("tests", "trackhub", package = "rtracklayer"))
uri(th)
genome(th)
length(th)
th$hg19
th[["hg19"]]
hub(th)
email(th)

\dontrun{
hub(th) <- "new_hub"
writeTrackHub(th)
}
}

\keyword{methods}
\keyword{classes}
