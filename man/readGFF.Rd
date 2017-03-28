\name{readGFF}

\alias{readGFF}
\alias{GFFcolnames}


\title{Reads a file in GFF format}

\description{
  Reads a file in GFF format and creates a data frame or
  \link[S4Vectors]{DataFrame} object from it. This is a low-level
  function that should not be called by user code.
}

\usage{
readGFF(filepath, version=0,
        columns=NULL, tags=NULL, filter=NULL, nrows=-1,
        raw_data=FALSE)

GFFcolnames(GFF1=FALSE)
}

\arguments{
  \item{filepath}{
    A single string containing the path or URL to the file to read.
    Alternatively can be a connection.
  }
  \item{version}{
    \code{readGFF} should do a pretty descent job at detecting the GFF
    version. Use this argument \emph{only} if it doesn't or if you want to
    force it to parse and import the file as if its 9-th column was in a
    different format than what it really is (e.g. specify \code{version=1}
    on a GTF or GFF3 file to interpret its 9-th column as the \code{"group"}
    column of a GFF1 file). Supported versions are 1, 2, and 3.
  }
  \item{columns}{
    The standard GFF columns to load. All of them are loaded by default.
  }
  \item{tags}{
    The tags to load. All of them are loaded by default.
  }
  \item{filter}{
  }
  \item{nrows}{
    \code{-1} or the maximum number of rows to read in (after filtering).
  }
  \item{raw_data}{
  }
  \item{GFF1}{
  }
}

\value{
  A DataFrame with columns corresponding to those in the GFF.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{import}} for importing a GFF file as a
          \link[GenomicRanges]{GRanges} object.

    \item \code{\link[GenomicRanges]{makeGRangesFromDataFrame}} in the
          \pkg{GenomicRanges} package for making a
          \link[GenomicRanges]{GRanges} object from a data frame or
          \link[S4Vectors]{DataFrame} object.

    \item \code{\link[GenomicFeatures]{makeTxDbFromGFF}} in the
          \pkg{GenomicFeatures} package for importing a GFF file as a
          \link[GenomicFeatures]{TxDb} object.

    \item The \link[S4Vectors]{DataFrame} class in the \pkg{S4Vectors}
          package.
  }
}

\examples{
## Standard GFF columns.
GFFcolnames()
GFFcolnames(GFF1=TRUE)  # "group" instead of "attributes"

tests_dir <- system.file("tests", package="rtracklayer")
test_gff3 <- file.path(tests_dir, "genes.gff3")

## Load everything.
df0 <- readGFF(test_gff3)
head(df0)

## Load some tags only (in addition to the standard GFF columns).
my_tags <- c("ID", "Parent", "Name", "Dbxref", "geneID")
df1 <- readGFF(test_gff3, tags=my_tags)
head(df1)

## Load no tags (in that case, the "attributes" standard column
## is loaded).
df2 <- readGFF(test_gff3, tags=character(0))
head(df2)

## Load some standard GFF columns only (in addition to all tags).
my_columns <- c("seqid", "start", "end", "strand", "type")
df3 <- readGFF(test_gff3, columns=my_columns)
df3
table(df3$seqid, df3$type)
makeGRangesFromDataFrame(df3, keep.extra.columns=TRUE)

## Combine use of 'columns' and 'tags' arguments.
readGFF(test_gff3, columns=my_columns, tags=c("ID", "Parent", "Name"))
readGFF(test_gff3, columns=my_columns, tags=character(0))

## Use the 'filter' argument to load only features of type "gene"
## or "mRNA" located on chr10.
my_filter <- list(type=c("gene", "mRNA"), seqid="chr10")
readGFF(test_gff3, filter=my_filter)
readGFF(test_gff3, columns=my_columns, tags=character(0), filter=my_filter)
}

\keyword{manip}
