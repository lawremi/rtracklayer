\name{liftOver}

\alias{liftOver}
\alias{liftOver,GenomicRanges,Chain-method}

\title{
  Lift intervals between genome builds
}
\description{
  A reimplementation of the UCSC liftover tool for lifting features from
  one genome build to another. In our preliminary tests, it is
  significantly faster than the command line tool. Like the UCSC tool, a
  chain file is required input.
}
\usage{
liftOver(x, chain, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The intervals to lift-over, usually
    a \code{\link[GenomicRanges:GRanges]{GRanges}}.
  }
  \item{chain}{
    A \code{\linkS4class{Chain}} object, usually imported with
    \code{\link{import.chain}}. 
  }
  \item{\dots}{
    Arguments for methods.
  }
}
\value{
  A \code{GRangesList} object. Each element contains the ranges mapped
  from the corresponding element in the input (may be one-to-many).
}
\references{
  \url{http://genome.ucsc.edu/cgi-bin/hgLiftOver}
}
\examples{
\dontrun{
chain <- import.chain("hg19ToHg18.over.chain")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx_hg18 <- liftOver(tx_hg19, chain)
}
}
\author{
  Michael Lawrence
}
