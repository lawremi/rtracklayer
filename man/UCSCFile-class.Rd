\name{UCSCFile-class}
\docType{class}

%% Classes:
\alias{class:UCSCFile}
\alias{UCSCFile-class}

%% Constructor:
\alias{UCSCFile}

%% Import:
\alias{import,UCSCFile,ANY,ANY-method}
\alias{import.ucsc}
\alias{import.ucsc,ANY-method}
\alias{import.ucsc,RTLFile-method}

%% Export:
\alias{export,ANY,UCSCFile,ANY-method}
\alias{export,RangedData,UCSCFile,ANY-method}
\alias{export,RangedDataList,UCSCFile,ANY-method}
\alias{export,UCSCData,UCSCFile,ANY-method}
\alias{export.ucsc}
\alias{export.ucsc,ANY,ANY-method}
\alias{export.ucsc,ANY,RTLFile-method}

\title{UCSCFile objects}

\description{
  These functions support the import and export of tracks emucscded
  within the UCSC track line metaformat, whereby multiple tracks may be
  concatenated within a single file, along with metadata mostly oriented
  towards visualization. Any \code{\linkS4class{UCSCData}} or
  \code{RangedDataList} object is
  automatically exported in this format, if the targeted format is known
  to be compatible. The BED and WIG import methods check for a track
  line, and delegate to these functions if one is found. Thus, calling
  this API directly is only necessary when importing embedded GFF
  (rare), or when one wants to create the track line during the export
  process.
}

\usage{
\S4method{import}{UCSCFile,ANY,ANY}(con, format, text,
                   subformat = "auto", drop = FALSE,
                   asRangedData = TRUE, genome = NA, ...)
import.ucsc(con, ...)

\S4method{export}{ANY,UCSCFile,ANY}(object, con, format, ...)
\S4method{export}{RangedData,UCSCFile,ANY}(object, con, format, ...)
\S4method{export}{RangedDataList,UCSCFile,ANY}(object, con, format,
                   append = FALSE, index = FALSE, ...)
\S4method{export}{UCSCData,UCSCFile,ANY}(object, con, format,
                   subformat = "auto", append = FALSE, index = FALSE, ...)
export.ucsc(object, con, ...)
}

\arguments{
  \item{con}{A path, URL, connection or \code{UCSCFile} object. For the
    functions ending in \code{.ucsc}, the file format is indicated by
    the function name. For the base \code{export} and \code{import}
    functions, \dQuote{ucsc} must be passed as the \code{format}
    argument.
  }
  \item{object}{The object to export, should be a \code{RangedData} or
    something coercible to a \code{RangedData}, like a
    \code{GRanges}. For exporting multiple tracks pass a
    \code{RangedDataList}, or something coercible to one, like a
    \code{GenomicRangesList}.
  }
  \item{format}{If not missing, should be \dQuote{ucsc}.}
  \item{text}{If \code{con} is missing, a character vector to use as the
    input
  }
  \item{subformat}{The file format to use for the actual features,
    between the track lines. Must be a text-based format that is
    compatible with track lines (most are). If an \code{RTLFile}
    subclass other than \code{UCSCFile} is passed as \code{con} to
    \code{import.ucsc} or \code{export.ucsc}, the subformat is assumed
    to be the corresponding format of \code{con}. Otherwise it defaults
    to \dQuote{auto}.  The following describes the logic of the
    \dQuote{auto} mode. For import, the subformat is taken as
    the \code{type} field in the track line. If none, the file
    extension is consulted.  For export, if \code{object} is a
    \code{UCSCData}, the subformat is taken as the \code{type}
    in its track line, if present. Otherwise, the subformat is chosen
    based on whether \code{object} contains a \dQuote{score} column. If
    there is a score, the target is either \code{BEDGraph} or
    \code{WIG}, depending on the structure of the ranges. Otherwise,
    \code{BED} is the target.
  }
  \item{genome}{The identifier of a genome, or \code{NA} if
    unknown. Typically, this is a UCSC identifier like \dQuote{hg19}. An
    attempt will be made to derive the \code{seqinfo} on the return
    value using either an installed BSgenome package or UCSC, if network
    access is available. This defaults to the \code{db} BED track line
    parameter, if any.
  }
  \item{asRangedData}{If \code{FALSE}, a \code{GenomicRangesList} is
    returned, instead of a \code{RangedDataList} of \code{UCSCData}
    objects, and the track line is then stored in the \code{metadata} of
    each \code{GRanges} in the list.
  }
  \item{drop}{If \code{TRUE}, and there is only one track in the file,
    return the track object directly, rather than embedding it in a list.
  }
  \item{append}{If \code{TRUE}, and \code{con} points to a file path,
    the data is appended to the file. Obviously, if \code{con} is a
    connection, the data is always appended.
  }
  \item{index}{If \code{TRUE}, automatically compress and index the
    output file with bgzf and tabix. Note that tabix indexing will
    sort the data by chromosome and start. Does not work when exporting
    a \code{RangedDataList} with multiple elements; tabix supports a
    single track in a file.
  }
  \item{...}{Should either specify track line parameters or arguments to
    pass down to the import and export routine for the subformat. 
  }
}

\value{
  A \code{RangedDataList} (or \code{GenomicRangesList} if the
  \code{asRangedData} is \code{FALSE}), unless \code{drop}
  is \code{TRUE} and there is only a single track in the file. In that
  case, the first and only object is extracted from the list and
  returned. The structure of that object depends on the format of the
  data. The \code{RangedDataList} contains \code{UCSCData} objects,
  whereas the \code{GenomicRangesList} contains \code{GRanges} objects
  with a \code{trackLine} element in their metadata.
}

\details{
  The UCSC track line permits the storage of multiple tracks in a single
  file by separating them with a so-called \dQuote{track line}, a line
  belonging with the word \dQuote{track} and containing
  various \code{key=value} pairs encoding metadata, most related to
  visualization. The standard fields in a track depend on the type of
  track being annotated. See \code{\linkS4class{TrackLine}} and its
  derivatives for how these lines are represented in R. The
  class \code{\linkS4class{UCSCData}} is an extension
  of \code{RangedData} with a formal slot for a \code{TrackLine}, and
  objects of this type are returned in a \code{RangedDataList} by the
  import functions described above, as long as \code{asRangedData}
  is \code{TRUE}. If \code{FALSE}, then each \code{GRanges} in the
  returned \code{GenomicRangesList} has the track line stored in its
  metadata, under the \code{trackLine} key. 

  For each track object to be exported, if the object is not a
  \code{UCSCData}, and there is no \code{trackLine} element in the
  metadata, then a new track line needs to be generated. This happens
  through the coercion of \code{object} to \code{UCSCData}. The track line
  is initialized to have the appropriate \code{type} parameter for the
  subformat, and the required \code{name} parameter is taken from the
  name of the track in the input list (if any). Otherwise, the default
  is simply \dQuote{R Track}. The \code{db} parameter (specific to BED
  track lines) is taken as \code{genome(object)} if not
  \code{NA}. Additional arguments passed to the export routines
  override parameters in the provided track line.

  If the subformat is either WIG or BEDGraph, and the features are
  stranded, a separate track will be output in the file for each
  strand. Neither of those formats encodes the strand and disallow
  overlapping features (which might occur upon destranding).
}

\section{UCSCFile objects}{
  The \code{UCSCFile} class extends \code{\linkS4class{RTLFile}} and is a
  formal represention of a resource in the UCSC format.
  To cast a path, URL or connection to a \code{UCSCFile}, pass it to
  the \code{UCSCFile} constructor.
}

\author{Michael Lawrence}

\references{
  \url{http://genome.ucsc.edu/goldenPath/help/customTrack.html}
}

\keyword{methods}
\keyword{classes}
