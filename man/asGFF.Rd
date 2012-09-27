\name{asGFF}
\alias{asGFF}
\alias{asGFF,GRangesList-method}

\title{Coerce to GFF structure}
\description{
  Coerce the structure of an object to one following GFF-like
  conventions, i.e., using the \code{Parent} GFF3 attribute to encode
  the hierarchical structure. This object is then suitable for export as GFF3.
}
\usage{
asGFF(x, ...)
\S4method{asGFF}{GRangesList}(x, parentType = "mRNA", childType = "exon")
}

\arguments{
  \item{x}{
    Generally, a tabular object to structure as GFF(3)
  }
  \item{parentType}{The value to store in the \code{type} column for the
    top-level (e.g., transcript) ranges.
  }
  \item{childType}{The value to store in the \code{type} column for the
    child (e.g., exon) ranges.
  }
  \item{\dots}{
    Arguments to pass to methods
  }
}

\value{
  For the \code{GRangesList} method:
  A \code{GRanges}, with the columns: \code{ID} (unique identifier),
  \code{Name} (from \code{names(x)}, and the names on each element of
  \code{x}, if any), \code{type} (as given by \code{parentType} and
  \code{childType}), and \code{Parent} (to relate each child range to
  its parent at the top-level).
}

\author{
  Michael Lawrence
}

\examples{
\dontrun{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exons <- exonsBy(TxDb_Hsapiens_UCSC_hg19_knownGene)
mcols(asGFF(exons))
}
}
