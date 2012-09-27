\name{QuickloadGenome-class}
\docType{class}

%% Classes:
\alias{class:QuickloadGenome}
\alias{QuickloadGenome-class}

%% Constructor:
\alias{QuickloadGenome}

%% Accessors:
\alias{seqinfo,QuickloadGenome-method}
\alias{seqinfo,DNAStringSet-method} % sneak it in here
\alias{seqinfo<-,QuickloadGenome-method}
\alias{genome,QuickloadGenome-method}
\alias{length,QuickloadGenome-method}
\alias{names,QuickloadGenome-method}
\alias{quickload}
\alias{mcols,QuickloadGenome-method}
\alias{releaseDate,QuickloadGenome-method}
\alias{organism,QuickloadGenome-method}
\alias{uri,QuickloadGenome-method}

%% Data Access:
\alias{track,QuickloadGenome-method}
\alias{track<-,QuickloadGenome,ANY-method}
\alias{track<-,QuickloadGenome,RangedData-method}
\alias{track<-,QuickloadGenome,RTLFile-method}
\alias{track<-,QuickloadGenome,RsamtoolsFile-method}
\alias{track<-,QuickloadGenome,character-method}
\alias{referenceSequence}
\alias{referenceSequence<-}
\alias{referenceSequence,QuickloadGenome-method}
\alias{referenceSequence<-,QuickloadGenome-method}

%% Show
\alias{show,QuickloadGenome-method}

\title{Quickload Genome Access}

\description{
  A Quickload data source is a collection of tracks and sequences,
  separated by genome. This class, \code{QuickloadGenome} provides
  direct access to the data for one particular genome.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{QuickloadGenome(quickload, genome, create = FALSE,
        seqinfo = seqinfo(genome),
        title = toString(genome))}: Constructs a
        new \code{QuickloadGenome} object, representing \code{genome} in
        the repository \code{quickload} (a URI string or a
        \code{\linkS4class{Quickload}} object).

        The \code{genome} argument can be an ID corresponding to a
        genome (potentially) in \code{quickload} or an installed
        \code{BSgenome} package. It can also be any instance of a class
        which has methods for \code{organism} and \code{releaseDate}. A
        good example is \code{\link[BSgenome:BSgenome-class]{BSgenome}}
        or any other derivative of
        \code{\link[BSgenome:GenomeDescription-class]{GenomeDescription}}.
        Those items are necessary for constructing the canonical
        Quickload genome string (G_Species_Month_Year).

        If \code{create} is \code{TRUE}, and the genome does not already
        exist, the genome will be created, using \code{seqinfo} for the
        sequence lengths and \code{title} for the display name of the
        genome in a UI. Creation only works if the repository is local
        and writeable. Reasonable defaults are used for \code{seqinfo}
        and \code{title} when the necessary methods are available (and
        they are for \code{BSgenome}).
    }
  }
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{Quickload}
  object.
  
  \describe{
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}: Gets or sets the
      \code{\link[GenomicRanges:Seqinfo-class]{Seqinfo}} object
      indicating the lengths of the sequences in the genome. No
      circularity information or genome identifier is stored.
    }
    \item{}{
      \code{quickload(x)}: Get the Quickload object that contains this
      genome.
    }
    \item{}{
      \code{uri(x)}: Get the uri pointing to the genome directory in the
      Quickload repository
    }
    \item{}{
      \code{genome(x)}: Get the name of the genome, e.g.
      \dQuote{H_sapiens_Feb_2009}.
    }
    \item{}{
      \code{releaseDate(x)}: Get the release portion of the genome name,
      e.g., \dQuote{Feb_2009}.
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
      \code{length(x)}: number of datasets
    }
    \item{}{
      \code{names(x), trackNames(x)}: names of the datasets
    }
    \item{}{
      \code{mcols(x)}: merged metadata on the datasets
    }
    \item{}{
      \code{track(x, name), x$name}: get the track called \code{name}
    }
    \item{}{
      \code{track(x, name, format = bestFileFormat(value), ...) <-
        value, x$name <- value}: store the track \code{value} under
        \code{name}. Note that track storing is only supported
        for local repositories, i.e., those with a \code{file://} URI
        scheme.

        Currently, supported \code{value} types include a
        \code{GenomicRanges}, \code{GRangesList}, or a file resource
        (copied to the repository).  The file resource may be
        represented as a path, URL, \code{\linkS4class{RTLFile}} or
        \code{\link[Rsamtools:RsamtoolsFile-class]{RsamtoolsFile}}. If
        not a file name, \code{value} is written in \code{format}. For
        generic interval data, this means a BigWig file (if there is a
        numeric \dQuote{score} column) or a BED file otherwise. An
        \code{RleList} (e.g., coverage) is output as BigWig. For
        \code{UCSCData} values, the format is chosen according to the
        type of track line. For \code{RsamtoolsFile} objects, the file
        and its index are copied.

        The arguments in \code{...} become attributes in the XML
        metadata. The \dQuote{description} attribute is standard and is
        a blurb for describing the track in a UI. For the rest, the
        interpretation is up to the client. IGB supports an ever-growing
        list; please see its documentation.
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
ql <- Quickload(file.path(tests_dir, "quickload"))
qlg <- QuickloadGenome(ql, "T_species_Oct_2011")
seqinfo(qlg)
organism(qlg)
releaseDate(qlg)
names(qlg)
mcols(qlg)
if (.Platform$OS.type != "windows") { # temporary
qlg$bedData
}

\dontrun{
## populating the test repository
ql <- Quickload(file.path(tests_dir, "quickload"), create = TRUE)
reference_seq <- import(file.path(tests_dir, "test.2bit"))
names(reference_seq) <- "test"
qlg <- QuickloadGenome(ql, "T_species_Oct_2011", create = TRUE,
                       seqinfo = seqinfo(reference_seq))
referenceSequence(qlg) <- reference_seq
test_bed <- import(file.path(tests_dir, "test.bed"))
names(test_bed) <- "test"
qlg$bedData <- test_bed
test_bedGraph <- import(file.path(tests_dir, "test.bedGraph"))
names(test_bedGraph) <- "test"
start(test_bedGraph) <- seq(1, 90, 10)
width(test_bedGraph) <- 10
track(qlg, "bedGraphData", format = "bw") <- test_bedGraph
}
}

\keyword{methods}
\keyword{classes}
