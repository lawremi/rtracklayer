\name{BEDFile-class}
\docType{class}

%% Classes:
\alias{class:BEDFile}
\alias{class:BED15File}
\alias{class:BEDGraphFile}
\alias{BEDFile-class}
\alias{BED15File-class}
\alias{BEDGraphFile-class}

%% Constructor:
\alias{BEDFile}
\alias{BED15File}
\alias{BEDGraphFile}

%% Import:
\alias{import,BEDFile,ANY,ANY-method}
\alias{import,BED15File,ANY,ANY-method}
\alias{import.bed}
\alias{import.bed,ANY-method}
\alias{import.bed15}
\alias{import.bed15,ANY-method}
\alias{import.bedGraph}
\alias{import.bedGraph,ANY-method}

%% Export:
\alias{export,ANY,BEDFile,ANY-method}
\alias{export,RangedData,BED15File,ANY-method}
\alias{export,RangedData,BEDFile,ANY-method}
\alias{export,RangedDataList,BEDFile,ANY-method}
\alias{export,UCSCData,BEDFile,ANY-method}
\alias{export.bed}
\alias{export.bed,ANY-method}
\alias{export.bed15}
\alias{export.bed15,ANY-method}
\alias{export.bedGraph}
\alias{export.bedGraph,ANY-method}

\title{BEDFile objects}

\description{
  These functions support the import and export of the UCSC BED
  format and its variants, including BEDGraph.
}

\usage{
\S4method{import}{BEDFile,ANY,ANY}(con, format, text, trackLine = TRUE,
                   genome = NA, asRangedData = TRUE, colnames = NULL,
                   which = NULL, seqinfo = NULL, extraCols = character())
import.bed(con, ...)
import.bed15(con, ...)
import.bedGraph(con,  ...)

\S4method{export}{ANY,BEDFile,ANY}(object, con, format, ...)
\S4method{export}{RangedData,BEDFile,ANY}(object, con, format,
                   append = FALSE, index = FALSE)
\S4method{export}{RangedDataList,BEDFile,ANY}(object, con, format, ...)
\S4method{export}{UCSCData,BEDFile,ANY}(object, con, format,
                   trackLine = TRUE, ...)
\S4method{export}{RangedData,BED15File,ANY}(object, con, format,
                   expNames = NULL, trackLine = NULL, ...)
export.bed(object, con, ...)
export.bed15(object, con, ...)
export.bedGraph(object, con, ...)
}

\arguments{
  \item{con}{A path, URL, connection or \code{BEDFile} object. For the
    functions ending in \code{.bed}, \code{.bedGraph} and \code{.bed15},
    the file format is indicated by the function name. For the base
    \code{export} and \code{import} functions, the format must be
    indicated another way. If \code{con} is a path, URL or connection,
    either the file extension or the \code{format} argument needs to be
    one of \dQuote{bed}, \dQuote{bed15} or \dQuote{bedGraph}. Compressed
    files (\dQuote{gz}, \dQuote{bz2} and \dQuote{xz}) are handled transparently.
  }
  \item{object}{The object to export, should be a \code{RangedData} or
    something coercible to a \code{RangedData}, like a
    \code{GRanges}. If the object has a method for \code{asBED}, it is
    called prior to coercion. This makes it possible to export a
    \code{GRangesList} or \code{TranscriptDb} in a way that preserves
    the hierarchical structure. For exporting multiple tracks, in the
    UCSC track line metaformat, pass a \code{RangedDataList}, or
    something coercible to one, like a \code{GenomicRangesList}.
  }
  \item{format}{If not missing, should be one of \dQuote{bed}, \dQuote{bed15}
    or \dQuote{bedGraph}.}
  \item{text}{If \code{con} is missing, a character vector to use as the
    input}
  \item{trackLine}{Whether to parse/output a UCSC track line. An
    imported track line will be stored in a \code{\linkS4class{TrackLine}}
    object, as part of the returned \code{\linkS4class{UCSCData}}.
  }
  \item{genome}{The identifier of a genome, or \code{NA} if
    unknown. Typically, this is a UCSC identifier like \dQuote{hg19}. An
    attempt will be made to derive the \code{seqinfo} on the return
    value using either an installed BSgenome package or UCSC, if network
    access is available.
  }
  \item{asRangedData}{If \code{FALSE}, a \code{GRanges} is returned,
    instead of a \code{RangedData}.
  }
  \item{colnames}{A character vector naming the columns to parse. These
    should name columns in the result, not those in the BED spec, so
    e.g. specify \dQuote{thick}, instead of \dQuote{thickStart}.
  }
  \item{which}{A range data structure like \code{RangesList}
    or \code{GRanges}. Only the intervals in the file overlapping the
    given ranges are returned. This is much more efficient when the
    file is indexed with the tabix utility.
  }
  \item{index}{If \code{TRUE}, automatically compress and index the
    output file with bgzf and tabix. Note that tabix indexing will
    sort the data by chromosome and start. Does not work when exporting
    a \code{RangedDataList} with multiple elements; tabix supports a
    single track in a file.
  }
  \item{seqinfo}{If not \code{NULL}, the \code{Seqinfo} object to set on
    the result. If the \code{genome} argument is not \code{NA}, it must
    agree with \code{genome(seqinfo)}.
  }
  \item{extraCols}{Names of extra columns to read from the BED file. As
    BED does not encode column names, these are assumed to be the last
    columns in the file. This enables parsing of the various BEDX+Y
    formats.
  }
  \item{append}{If \code{TRUE}, and \code{con} points to a file path,
    the data is appended to the file. Obviously, if \code{con} is a
    connection, the data is always appended.
  }
  \item{expNames}{character vector of column names in \code{object} to
    export as sample columns in the BED15 file.}
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{BEDFile} method on
    \code{import}. For export, the \code{RangedData}, \code{BEDFile}
    method on \code{export} is the sink. When \code{trackLine} is
    \code{TRUE} or the target format is BED15, the arguments are passed
    through \code{export.ucsc}, so track line parameters are supported.
  }
}

\value{
  A \code{RangedData} (or \code{GRanges} if the \code{asRangedData} is
  \code{FALSE}), with the columns described in the details.
}

\details{
  The BED format is a tab-separated table of intervals, with annotations
  like name, score and even sub-intervals for representing alignments
  and gene models. Official (UCSC) child formats currently include BED15
  (adding a number matrix for e.g. expression data across multiple
  samples) and BEDGraph (a compressed means of storing a single score
  variable, e.g. coverage; overlapping features are not allowed). Many
  tools and organizations have extended the BED format with additional
  columns for particular use cases. These are not yet supported by
  rtracklayer, but a mechanism will be added soon. The advantage of BED
  is its balance between simplicity and expressiveness. It is also
  relatively scalable, because only the first three columns (chrom,
  start and end) are required. Thus, BED is best suited for representing
  simple features. For specialized cases, one is usually better off with
  another format. For example, genome-scale vectors belong in
  \link[=BigWigFile]{BigWig}, alignments from high-throughput sequencing
  belong in \link[Rsamtools:BamFile]{BAM}, and gene models are more
  richly expressed in \link[=GFFFile]{GFF}.

  The following is the mapping of BED elements to a \code{RangedData}
  object. NA values are allowed only where indicated. These appear as
  a \dQuote{.} in the file. Only the first three columns (chrom, start
  and strand) are required. The other columns can only be included if
  all previous columns (to the left) are included. Upon export, default
  values are used to automatically pad the table, if necessary.
  \describe{
    \item{chrom, start, end}{the \code{ranges} component.}
    \item{name}{character vector (NA's allowed) in the \code{name}
      column; defaults to NA on export. 
    }
    \item{score}{numeric vector in the \code{score}
      column, accessible via the \code{score} accessor. Defaults to 0
      on export. This is the only column present in BEDGraph (besides
      chrom, start and end), and it is required.
    }
    \item{strand}{strand factor (NA's allowed) in the \code{strand}
      column, accessible via the \code{strand} accessor; defaults to NA
      on export.
    }
    \item{thickStart, thickEnd}{\code{Ranges} object in a
      column named \code{thick}; defaults to the ranges of the feature
      on export.
    }
    \item{itemRgb}{character vector of hex color codes, as returned
      by \code{\link{col2rgb}}, in the \code{itemRgb} column; default is
      NA on export, which translates to black.
    }
    \item{blockSizes, blockStarts, blockCounts}{\code{RangesList} object
      in a column named \code{blocks}; defaults to empty upon BED15 export.
    }
  }

  These columns are present only in BED15:
  \describe{
    \item{expCount, expIds, expScores}{A column for each unique element
      in \code{expIds}, containing the corresponding values from
      \code{expScores}. When a value is not present for a feature, NA is
      substituted. NA values become -10000 in the file.
    }
  }
}

\section{BEDFile objects}{
  The \code{BEDFile} class extends \code{\linkS4class{RTLFile}} and is a
  formal represention of a resource in the BED format.
  To cast a path, URL or connection to a \code{BEDFile}, pass it to
  the \code{BEDFile} constructor. Classes and constructors also exist
  for the subclasses \code{BED15File} and \code{BEDGraphFile}.
}

\author{Michael Lawrence}

\references{
  \url{http://genome.ucsc.edu/goldenPath/help/customTrack.html}
}

\examples{
  test_path <- system.file("tests", package = "rtracklayer")
  test_bed <- file.path(test_path, "test.bed")

  test <- import(test_bed)
  test
  import.bed(test_bed)

  test_bed_file <- BEDFile(test_bed)
  import(test_bed_file)
  
  test_bed_con <- file(test_bed)
  import(test_bed_con, format = "bed")
  close(test_bed_con)
  
  import(test_bed, trackLine = FALSE)
  import(test_bed, asRangedData = FALSE)
  import(test_bed, genome = "hg19")
  import(test_bed, colnames = c("name", "strand", "thick"))
  
  which <- RangesList(chr7 = ranges(test)[[1]][1:2])
  import(test_bed, which = which)

\dontrun{
  test_bed_out <- file.path(tempdir(), "test.bed")
  export(test, test_bed_out)
  
  test_bed_out_file <- BEDFile(test_bed_out)
  export(test, test_bed_out_file)

  export(test, test_bed_out, name = "Alternative name")
  
  test_bed_gz <- paste(test_bed_out, ".gz", sep = "")
  export(test, test_bed_gz)
  
  export(test, test_bed_out, index = TRUE)
  export(test, test_bed_out, index = TRUE, trackLine = FALSE)
  
  bed_text <- export(test, format = "bed")
  test <- import(format = "bed", text = bed_text)
}
}

\keyword{methods}
\keyword{classes}
