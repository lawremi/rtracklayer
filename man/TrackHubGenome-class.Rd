\name{QuickloadGenome-class}
\docType{class}

%% Classes:
\alias{class:TrackHubGenome}
\alias{TrackHubGenome-class}

%% Constructor:
\alias{TrackHubGenome}

%% Accessors:
\alias{genome,TrackHubGenome-method}
\alias{length,TrackHubGenome-method}
\alias{names,TrackHubGenome-method}
\alias{trackhub}
\alias{uri,TrackHubGenome-method}

%% Data Access:
\alias{referenceSequence,TrackHubGenome-method}
\alias{referenceSequence<-,TrackHubGenome-method}

%% Show
\alias{show,TrackHubGenome-method}

\title{TrackHub Genome Access}

\description{
  A TrackHub data source is a collection of tracks and sequences,
  separated by genome. This class, \code{TrackHubGenome} provides
  direct access to the data for one particular genome.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{TrackHubGenome(trackhub, genome, create = FALSE,
        genomeRecord = Genome()}: Constructs a
        new \code{TrackHubGenome} object, representing \code{genome} in
        the repository \code{trackhub} (a URI string or a
        \code{\linkS4class{TrackHub}} object).

        The \code{genome} argument can be an ID corresponding to a
        genome (potentially) in \code{trackhub} or an installed
        \code{BSgenome} package.

        If \code{create} is \code{TRUE}, and the genome does not already
        exist, the genome will be created, using \code{Genome}.
        \code{Genome} takes names arguemnts of
        Creation only works if the repository is local and
        writeable.
    }
  }
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} and \code{object} represent a
  \code{TrackHubGenome} object.

  \describe{
    \item{}{
      \code{trackhub(x)}: Get the TracckHub object that contains this
      genome.
    }
    \item{}{
      \code{uri(x)}: Get the uri pointing to the genome directory in the
      TrackHub repository.
    }
    \item{}{
      \code{genome(x)}: Get the name of the genome, e.g.
      \dQuote{hg19}.
    }
    \item{}{
      \code{organism(x)}: Get the organism portion of the genome name,
      e.g., \dQuote{H sapiens}.
    }
  }
}

\section{Data Access}{
  \describe{
    \item{}{
      \code{referenceSequence(x)}: Get the reference sequence, as a
      \code{DNAStringSet}.
    }
    \item{}{
      \code{referenceSequence(x) <- value}: Set the reference sequence, as a
      \code{DNAStringSet}. It is written as a 2bit file. This only works
      on local repositories.
    }
  }
}

\author{Michael Lawrence}

\examples{
tests_dir <- system.file("tests", package = "rtracklayer")
th <- TrackHub(file.path(tests_dir, "trackhub"))
thg <- TrackHubGenome(th, "hg19")
length(thg)
organism(thg)

\dontrun{
th <- TrackHub(file.path(tests_dir, "trackhub"), create = TRUE)
genomeFile(th) <- "genomes.txt"
hg38 <- Genome(genome = "hg38")
thg <- TrackHubGenome(th, "hg38", create=TRUE, hg38)
}
}

\keyword{methods}
\keyword{classes}
