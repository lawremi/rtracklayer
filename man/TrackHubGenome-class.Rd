\name{TrackHubGenome-class}
\docType{class}

%% Classes:
\alias{class:TrackHubGenome}
\alias{TrackHubGenome-class}

%% Constructor:
\alias{TrackHubGenome}

%% Accessors:
\alias{uri,TrackHubGenome-method}
\alias{genome,TrackHubGenome-method}
\alias{length,TrackHubGenome-method}
\alias{names,TrackHubGenome-method}
\alias{getTracks,TrackHubGenome-method}
\alias{trackhub}
\alias{organism,TrackHubGenome-method}
\alias{trackField,TrackHubGenome-method}
\alias{trackField<-,TrackHubGenome-method}
\alias{writeTrackHub,TrackHubGenome-method}

%% Data Access:
\alias{track,TrackHubGenome-method}
\alias{track<-,TrackHubGenome,ANY-method}
\alias{track<-,TrackHubGenome,RTLFile-method}
\alias{track<-,TrackHubGenome,RsamtoolsFile-method}
\alias{track<-,TrackHubGenome,character-method}
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
      \code{TrackHubGenome(trackhub, genome, create = FALSE}: Constructs a
        new \code{TrackHubGenome} object, representing \code{genome} in
        the repository \code{trackhub} (a URI string or a
        \code{\linkS4class{TrackHub}} object).

        The \code{genome} argument can be an ID corresponding to a
        genome (potentially) in \code{trackhub} or an installed
        \code{BSgenome} package.

        If \code{create} is \code{TRUE}, and the trackDb file does not
        already exist, it will be created.
        Creation only works if the repository is local and
        writeable.
    }
  }
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represent a \code{TrackHubGenome} object.

  \describe{
    \item{}{
      \code{uri(x)}: Get the uri pointing to the genome directory in the
      TrackHub repository.
    }
    \item{}{
      \code{genome(x)}: Get the name of the genome, e.g.
      \dQuote{hg19}.
    }
    \item{}{
      \code{length(x)}: number of tracks
    }
    \item{}{
      \code{names(x), trackNames(x)}: names of the tracks
    }
    \item{}{
      \code{getTracks(x)}: Get the \code{List} of \code{Track} from the tracks
    }
    \item{}{
      \code{trackhub(x)}: Get the TrackHub object that contains this
      genome.
    }
    \item{}{
      \code{organism(x)}: Get the organism name for this genome,
      e.g., \dQuote{H sapiens}.
    }
    \item{}{
      \code{trackField(x, name, field)}: Get the \code{value} of \code{field} for \code{name} track.
    }
    \item{}{
      \code{trackField(x, name, field) <- value}: Store the \code{field} and \code{value} for \code{name} track.
    }
    \item{}{\code{writeTrackHub(x)}:
      Write tracks from memory representation to the trackDb file.
    }
  }
}

\section{Data Access}{
  \describe{
    \item{}{
      \code{track(x, name), x$name}: get the track called \code{name}
    }
    \item{}{
      \code{track(x, name, format = bestFileFormat(value)) <-
        value, x$name <- value}: store the track \code{value} under
        \code{name}. Note that track storing is only supported
        for local repositories, i.e., those with a \code{file://} URI
        scheme.

        Currently, supported \code{value} types include a
        \code{GenomicRanges}, \code{GRangesList}, or a file resource
        (copied to the repository).  The file resource may be
        represented as a path, URL, \code{\linkS4class{BiocFile}} or
        \code{\link[Rsamtools:RsamtoolsFile-class]{RsamtoolsFile}}. If
        not a file name, \code{value} is written in \code{format}. For
        generic interval data, this means a BigWig file (if there is a
        numeric \dQuote{score} column) or a BED file otherwise. An
        \code{RleList} (e.g., coverage) is output as BigWig. For
        \code{UCSCData} values, the format is chosen according to the
        type of track line. For \code{RsamtoolsFile} objects, the file
        and its index are copied.
    }
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
names(thg)

\dontrun{
th <- TrackHub(file.path(tests_dir, "trackhub"), create = TRUE)
genomesFile(th) <- "genomes.txt"
genomeInfo(th) <- Genome(genome = "hg38", trackDb = "hg38/trackDb.txt")
genomeField(th, "hg38", "twoBitPath") <- "hg38/seq.2bit"
writeTrackHub(th)
thg <- TrackHubGenome(th, "hg38", create = TRUE)
seq <- import(file.path(tests_dir, "test.2bit"))
referenceSequence(thg) <- seq
track(thg, "PeaksData") <- paste0(tests_dir, "/test.bigWig")
trackField(thg, "wgEncodeUWDukeDnaseGM12878FdrPeaks", "bigDataUrl") <- "hg38/wgEncodeCshlShortRnaSeq.bigWig"
trackField(thg, "wgEncodeUWDukeDnaseGM12878FdrPeaks", "color") <- "8,104,172"
writeTrackHub(thg)
}
}

\keyword{methods}
\keyword{classes}
