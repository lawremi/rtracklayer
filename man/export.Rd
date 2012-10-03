\name{io}

\alias{export}
\alias{export,ANY,connection,character-method}
\alias{export,ANY,connection,missing-method}
\alias{export,ANY,CompressedFile,missing-method}
\alias{export,ANY,missing,character-method}
\alias{export,ANY,character,missing-method}
\alias{export,ANY,character,character-method}

\alias{import}
\alias{import,connection,character,ANY-method}
\alias{import,connection,missing,ANY-method}
\alias{import,CompressedFile,missing,ANY-method}
\alias{import,character,missing,ANY-method}
\alias{import,character,character,ANY-method}
\alias{import,missing,ANY,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Import and export  }
\description{
  The functions \code{import} and \code{export} load and save
  objects from and to particular file formats. The rtracklayer package
  implements support for a number of annotation and sequence formats.
}
\usage{
export(object, con, format, ...)
import(con, format, text, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The object to export. }
  \item{con}{ The connection from which data is loaded or to which data
    is saved. If this is a character vector, it is assumed to be a
    filename and a corresponding file connection is created and then
    closed after exporting the object. If a \code{\linkS4class{RTLFile}}
    derivative, the data is loaded from or saved to the underlying
    resource. If missing, the function will return the output as a
    character vector, rather than writing to a connection.
  }
  \item{format}{ The format of the output. If missing and \code{con} is
    a filename, the format is derived from the file extension. This
    argument is unnecessary when \code{con} is a derivative of
    \code{\linkS4class{RTLFile}}.
  }
  \item{text}{ If \code{con} is missing, this can be a character vector
    directly providing the string data to import. }
  \item{\dots}{ Parameters to pass to the format-specific method.
  }
}

\details{
  The rtracklayer package supports a number of file formats for
  representing annotated genomic intervals. These are
  each represented as a subclass of \code{\linkS4class{RTLFile}}. Below,
  we list the major supported formats, with some advice for when a
  particular file format is appropriate:
  
  \describe{
    \item{\acronym{\link[=GFFFile]{GFF}}}{The General Feature Format is
      meant to represent any set of genomic features, with
      application-specific columns represented as
      \dQuote{attributes}. There are three principal versions (1, 2, and
      3). This is a good format for interoperating with other genomic
      tools and is the most flexible format, in that a feature may have
      any number of attributes (in version 2 and above). Version 3
      (GFF3) is the preferred version. Its specification lays out
      conventions for representing various types of data, including gene
      models, for which it is the format of choice. For variants,
      rtracklayer has rudimentary support for an extention of GFF3
      called GVF. UCSC supports GFF1, but it needs to be encapsulated in
      the UCSC metaformat, i.e. \code{export.ucsc(subformat =
      "gff1")}. The BED format is typically preferred over GFF for
      interaction with UCSC. GFF files can be indexed with the tabix
      utility for fast range-based queries via rtracklayer and
      Rsamtools.
    }
    
    \item{\acronym{\link[=BEDFile]{BED}}}{The Browser Extended Display
      format is for displaying qualitative tracks in a genome browser,
      in particular UCSC. It finds a good balance between simplicity and
      expressiveness. It is much simpler than GFF and yet can still
      represent multi-exon gene structures. It is somewhat limited by
      its lack of the attribute support of GFF. To circumvent this, many
      tools and organizations have extended BED with additional
      columns. These are not officially valid BED files, and as such
      rtracklayer does not yet support them (this will be addressed
      soon). The rtracklayer package does support two official
      extensions of BED: Bed15 and bedGraph, see below.  BED files can
      be indexed with the tabix utility for fast range-based queries via
      rtracklayer and Rsamtools.
    }
 
    \item{\acronym{\link[=BED15File]{Bed15}}}{An extension of BED with 15
      columns, Bed15 is meant to represent data from microarray
      experiments. Multiple samples/columns are supported, and the data
      is displayed in UCSC as a compact heatmap. Few other tools support
      this format. With 15 columns per feature, this format is probably
      too verbose for e.g. ChIP-seq coverage (use multiple BigWig tracks
      instead).}
    
    \item{\acronym{\link[=BEDGraphFile]{bedGraph}}}{A variant of BED that
      represents a score column more compactly than \acronym{BED} and
      especially \acronym{Bed15}, although only one sample is
      supported. The data is displayed in UCSC as a bar or line
      graph. For large data (the typical case), \code{BigWig} is
      preferred.
    }
    
    \item{\acronym{\link[=WIGFile]{WIG}}}{The Wiggle format is meant for
      storing dense numerical data, such as window-based GC and
      conservation scores. The data is displayed in UCSC as a bar or
      line graph. The WIG format only works for intervals with a uniform
      width. For non-uniform widths, consider \code{bedGraph}. For large
      data, consider \code{BigWig}.
    }

    \item{\acronym{\link[=BigWigFile]{BigWig}}}{The BigWig format is a
      binary version of both \code{bedGraph} and \code{WIG} (which are
      now somewhat obsolete). A BigWig file contains a spatial index for
      fast range-based queries and also embeds summary statistics of the
      scores at several zoom levels. Thus, it is ideal for visualization
      of and parallel computing on genome-scale vectors, like the
      coverage from a high-throughput sequencing experiment.
    }  
  }
  
  In summary, for the typical use case of combining gene models with
  experimental data, \acronym{GFF} is preferred for gene models and
  \code{BigWig} is preferred for quantitative score vectors. Note that
  the Rsamtools package provides support for the
  \code{\link[Rsamtools:BamFile]{BAM}} file format (for representing
  read alignments), among others. Based on this, the rtracklayer package
  provides an \code{export} method for writing \code{GappedAlignments}
  and \code{GappedReads} objects as \code{BAM}. For variants, consider
  VCF, supported by the VariantAnnotation package.

  There is also support for reading and writing biological sequences,
  including the UCSC \code{\link[=TwoBitFile]{TwoBit}} format for
  compactly storing a genome sequence along with a mask. The files are
  binary, so they are efficiently queried for particular ranges. A
  similar format is \code{\link[Rsamtools:FaFile]{FA}}, supported by
  Rsamtools.
}

\value{
  If \code{con} is missing, a character vector containing the string
  output. Otherwise, nothing is returned.
}
\author{ Michael Lawrence }

\examples{
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
  \dontrun{export(track, "my.gff", version = "3")}
  ## equivalently,
  \dontrun{export(track, "my.gff3")}
  ## or
  \dontrun{
  con <- file("my.gff3")
  export(track, con, "gff3")
  close(con)
  }
  ## or as a string
  export(track, format = "gff3")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
