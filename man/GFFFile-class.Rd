\name{GFFFile-class}
\docType{class}

%% Classes:
\alias{class:GFFFile}
\alias{class:GFF1File}
\alias{class:GFF2File}
\alias{class:GFF3File}
\alias{class:GVFFile}
\alias{class:GTFFile}
\alias{GFFFile-class}
\alias{GFF1File-class}
\alias{GFF2File-class}
\alias{GFF3File-class}
\alias{GVFFile-class}
\alias{GTFFile-class}

%% Constructor:
\alias{GFFFile}
\alias{GFF1File}
\alias{GFF2File}
\alias{GFF3File}
\alias{GVFFile}
\alias{GTFFile}

%% Import:
\alias{import,GFFFile,ANY,ANY-method}
\alias{import.gff}
\alias{import.gff1}
\alias{import.gff2}
\alias{import.gff3}
\alias{import.gff,ANY-method}
\alias{import.gff1,ANY-method}
\alias{import.gff2,ANY-method}
\alias{import.gff3,ANY-method}

%% Export:
\alias{export,ANY,GFFFile,ANY-method}
\alias{export,RangedData,GFFFile,ANY-method}
\alias{export,RangedDataList,GFFFile,ANY-method}
\alias{export.gff}
\alias{export.gff,ANY-method}
\alias{export.gff1}
\alias{export.gff1,ANY-method}
\alias{export.gff2}
\alias{export.gff2,ANY-method}
\alias{export.gff3}
\alias{export.gff3,ANY-method}

% Other:
\alias{genome,GFFFile-method}

\title{GFFFile objects}

\description{
  These functions support the import and export of the GFF format, of
  which there are three versions and several flavors.
}

\usage{
\S4method{import}{GFFFile,ANY,ANY}(con, format, text,
           version = c("", "1", "2", "3"),
           genome = NA, asRangedData = TRUE, colnames = NULL, which = NULL,
           feature.type = NULL)
import.gff(con, ...)
import.gff1(con, ...)
import.gff2(con, ...)
import.gff3(con, ...)

\S4method{export}{ANY,GFFFile,ANY}(object, con, format, ...)
\S4method{export}{RangedData,GFFFile,ANY}(object, con, format,
                   version = c("1", "2", "3"),
                   source = "rtracklayer", append = FALSE, index = FALSE)
\S4method{export}{RangedDataList,GFFFile,ANY}(object, con, format, ...)
export.gff(object, con, ...)
export.gff1(object, con, ...)
export.gff2(object, con, ...)
export.gff3(object, con, ...)
}

\arguments{
  \item{con}{A path, URL, connection or \code{GFFFile} object. For the
    functions ending in \code{.gff}, \code{.gff1}, etc, the file format
    is indicated by the function name. For the base \code{export} and
    \code{import} functions, the format must be indicated another
    way. If \code{con} is a path, URL or connection, either the file
    extension or the \code{format} argument needs to be one of
    \dQuote{gff}, \dQuote{gff1} \dQuote{gff2}, \dQuote{gff3},
    \dQuote{gvf}, or \dQuote{gtf}. Compressed files (\dQuote{gz},
    \dQuote{bz2} and \dQuote{xz}) are handled transparently.
  }
  \item{object}{The object to export, should be a \code{RangedData} or
    something coercible to a \code{RangedData}, like a
    \code{GRanges}. If the object has a method for \code{asGFF}, it is
    called prior to coercion. This makes it possible to export a
    \code{GRangesList} or \code{TranscriptDb} in a way that preserves
    the hierarchical structure. For exporting multiple tracks, in the UCSC
    track line metaformat, pass a \code{RangedDataList}, or something
    coercible to one, like a \code{GenomicRangesList}.
  }
  \item{format}{If not missing, should be one of \dQuote{gff}, \dQuote{gff1}
    \dQuote{gff2}, \dQuote{gff3}, \dQuote{gvf}, or \dQuote{gtf}.
  }
  \item{version}{If the format is given as \dQuote{gff}, i.e., it does
    not specify a version, then this should indicate the GFF version as
    one of \dQuote{} (for import only, from the \code{gff-version}
    directive in the file or \dQuote{1} if none), \dQuote{1}, \dQuote{2}
    or \dQuote{3}.
  }
  \item{text}{If \code{con} is missing, a character vector to use as the
    input.
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
    should name either fixed fields, like \code{source} or
    \code{strand}, or, for GFF2 and GFF3, any attribute.
  }
  \item{which}{A range data structure like \code{RangesList}
    or \code{GRanges}. Only the intervals in the file overlapping the
    given ranges are returned. This is much more efficient when the
    file is indexed with the tabix utility.
  }
  \item{feature.type}{\code{NULL} (the default) or a character vector of
    valid feature types. If not \code{NULL}, then only the features of the
    specified type(s) are imported.
  }
  \item{source}{The value for the source column in GFF. This is
    typically the name of the package or algorithm that generated the
    feature.
  }
  \item{index}{If \code{TRUE}, automatically compress and index the
    output file with bgzf and tabix. Note that tabix indexing will
    sort the data by chromosome and start. Does not work when exporting
    a \code{RangedDataList} with multiple elements; tabix supports a
    single track in a file.
  }
  \item{append}{If \code{TRUE}, and \code{con} points to a file path,
    the data is appended to the file. Obviously, if \code{con} is a
    connection, the data is always appended.
  }
  \item{...}{Arguments to pass down to methods to other methods. For
    import, the flow eventually reaches the \code{GFFFile} method on
    \code{import}. For export, the \code{RangedData}, \code{GFFFile}
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
  The Generic Feature Format (GFF) format is a tab-separated table of
  intervals. There are three different versions of GFF, and they all
  have the same number of columns. In GFF1, the last column is a
  grouping factor, whereas in the later versions the last column holds
  application-specific attributes, with some conventions defined for
  those commonly used. This attribute support facilitates specifying
  extensions to the format. These include GTF (Gene Transfer Format, an
  extension of GFF2) and GVF (Genome Variation Format, an extension of
  GFF3).  The rtracklayer package recognizes the \dQuote{gtf} and
  \dQuote{gvf} extensions and parses the extra attributes into columns
  of the result; however, it does not perform any extension-specific
  processing. Both GFF1 and GFF2 have been proclaimed obsolete; however,
  the UCSC Genome Browser only supports GFF1 (and GTF), and GFF2 is still
  in broad use.

  GFF is distinguished from the simpler BED format by its flexible
  attribute support and its hierarchical structure, as specified by the
  \code{group} column in GFF1 (only one level of grouping) and the
  \code{Parent} attribute in GFF3. GFF2 does not specify a convention
  for representing hierarchies, although its GTF extension provides this
  for gene structures. The combination of support for hierarchical data
  and arbitrary descriptive attributes makes GFF(3) the preferred format
  for representing gene models.

  Although GFF features a \code{score} column, large quantitative data
  belong in a format like \link[=BigWigFile]{BigWig} and alignments from
  high-throughput experiments belong in
  \link[Rsamtools:BamFile]{BAM}. For variants, the VCF format (supported
  by the VariantAnnotation package) seems to be more widely adopted than
  the GVF extension.

  A note on the UCSC track line metaformat: track lines are a means for
  passing hints to visualization tools like the UCSC Genome Browser and
  the Integrated Genome Browser (IGB), and they allow multiple tracks to
  be concatenated in the same file. Since GFF is not a UCSC format, it
  is not common to annotate GFF data with track lines, but rtracklayer
  still supports it. To export or import GFF data in the track line
  format, call \code{\link{export.ucsc}} or \code{\link{import.ucsc}}.
  
  The following is the mapping of GFF elements to a \code{RangedData}
  object. NA values are allowed only where indicated. These appear as
  a \dQuote{.} in the file. GFF requires that all columns are included,
  so \code{export} generates defaults for missing columns.

  \describe{
    \item{seqid, start, end}{the \code{ranges} component.}
    \item{source}{character vector in the \code{source}
      column; defaults to \dQuote{rtracklayer} on export.
    }
    \item{type}{character vector in the \code{type} column; defaults
      to \dQuote{sequence_feature} in the output, i.e., SO:0000110.
    }
    \item{score}{numeric vector (NA's allowed) in the \code{score}
      column, accessible via the \code{score} accessor; defaults
      to \code{NA} upon export.
    }
    \item{strand}{strand factor (NA's allowed) in the \code{strand}
      column, accessible via the \code{strand} accessor; defaults
      to \code{NA} upon export.
    }
    \item{phase}{integer vector, either 0, 1 or 2 (NA's allowed);
      defaults to \code{NA} upon export.
    }
    \item{group}{a factor (GFF1 only); defaults to the \code{seqid}
      (e.g., chromosome) on export.
    }
  }

  In GFF versions 2 and 3, attributes map to arbitrary columns in the
  result. In GFF3, some attributes (\code{Parent}, \code{Alias},
  \code{Note}, \code{DBxref} and \code{Ontology_term}) can have
  multiple, comma-separated values; these columns are thus always
  \code{CharacterList} objects.
}

\section{GFFFile objects}{
  The \code{GFFFile} class extends \code{\linkS4class{RTLFile}} and is a
  formal represention of a resource in the GFF format.
  To cast a path, URL or connection to a \code{GFFFile}, pass it to
  the \code{GFFFile} constructor. The \code{GFF1File}, \code{GFF2File},
  \code{GFF3File}, \code{GVFFile} and \code{GTFFile} classes all extend
  \code{GFFFile} and indicate a particular version of the format.

  It has the following utility methods:
  \describe{
    \item{}{\code{genome}: Gets the genome identifier from
      the \dQuote{genome-build} header directive.
    }
  }
}

\author{Michael Lawrence}

\references{
  \describe{
    \item{GFF1, GFF2}{
      \url{http://www.sanger.ac.uk/resources/software/gff/spec.html}
    }
    \item{GFF3}{\url{http://www.sequenceontology.org/gff3.shtml}}
    \item{GVF}{\url{http://www.sequenceontology.org/resources/gvf.html}}
    \item{GTF}{\url{http://mblab.wustl.edu/GTF22.html}}
  }
}

\examples{
  test_path <- system.file("tests", package = "rtracklayer")
  test_gff3 <- file.path(test_path, "genes.gff3")

  ## basic import
  test <- import(test_gff3)
  test

  ## import.gff functions
  import.gff(test_gff3)
  import.gff3(test_gff3)

  ## GFFFile derivatives
  test_gff_file <- GFF3File(test_gff3)
  import(test_gff_file)
  test_gff_file <- GFFFile(test_gff3)
  import(test_gff_file)
  test_gff_file <- GFFFile(test_gff3, version = "3")
  import(test_gff_file)

  ## from connection
  test_gff_con <- file(test_gff3)
  test <- import(test_gff_con, format = "gff")
  close(test_gff_con)

  ## various arguments
  import(test_gff3, genome = "hg19")
  import(test_gff3, asRangedData = FALSE)
  import(test_gff3, colnames = character())
  import(test_gff3, colnames = c("strand", "geneName"))

  ## 'which'
  which <- RangesList(chr10 = IRanges(90000, 93000))
  import(test_gff3, which = which)

\dontrun{
  ## 'append'
  test_gff3_out <- file.path(tempdir(), "genes.gff3")
  export(test["chr10"], test_gff3_out)
  export(test["chr12"], test_gff3_out, append = TRUE)
  import(test_gff3_out)
  
  ## 'index'
  export(test, test_gff3_out, index = TRUE)
  test_bed_gz <- paste(test_gff3_out, ".gz", sep = "")
  import(test_bed_gz, which = which)
  
  ## RangedDataList  
  rdl <-
    RangedDataList(new("UCSCData", test[1],
                       trackLine = new("BasicTrackLine", name = "chr10")),
                   new("UCSCData", test[2],
                       trackLine = new("BasicTrackLine", name = "chr12")))
  names(rdl) <- names(test)
  export(rdl, test_gff3_out)
  import.ucsc(test_gff3_out)
}
}

\keyword{methods}
\keyword{classes}
