\name{asBED}
\alias{asBED}
\alias{asBED,GRangesList-method}

\title{Coerce to BED structure}
\description{
  Coerce the structure of an object to one following BED-like
  conventions, i.e., with columns for blocks and thick regions.
}
\usage{
asBED(x, ...)
\S4method{asBED}{GRangesList}(x)
}

\arguments{
  \item{x}{
    Generally, a tabular object to structure as BED
  }
  \item{\dots}{
    Arguments to pass to methods
  }
}

\details{
  The exact behavior depends on the class of \code{object}.
  \describe{
    \item{\code{GRangesList}}{This treats \code{object} as if it were a
      list of transcripts, i.e., each element contains the exons of a
      transcript. The \code{blockStarts} and
      \code{blockSizes} columns are derived from the ranges in each
      element. Also, add \code{name} column from \code{names(object)}.
    }
  }
}

\value{
  A \code{GRanges}, with the columns \code{name},
  \code{blockStarts} and \code{blockSizes} added.
}

\author{
  Michael Lawrence
}

\examples{
\dontrun{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exons <- exonsBy(TxDb_Hsapiens_UCSC_hg19_knownGene)
mcols(asBED(exons))
}
}
