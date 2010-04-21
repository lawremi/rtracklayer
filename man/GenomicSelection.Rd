\name{GenomicSelection}
\alias{GenomicSelection}
\title{
Genomic data selection
}
\description{
  Convenience constructor of a \link[IRanges]{RangedSelection} object
  for selecting a data on a per-chromosome basis for a given genome. 
}
\usage{
GenomicSelection(genome, chrom = NULL, colnames = character(0))
}
\arguments{
  \item{genome}{
    A string identifying a genome. Should match the end of a BSgenome
    package name, e.g. "hg19".
  }
  \item{chrom}{
    Character vector naming chromosomes to select.
  }
  \item{colnames}{
    The column names to select from the dataset.
  }
}
\value{
  A \link[IRanges]{RangedSelection} object, selecting entire chromosomes
}
\author{
  Michael Lawrence
}

\seealso{
  \code{\link[IRanges]{RangedSelection}}, \code{\link{BigWigSelection}}
}
\examples{
# every chromosome from hg19
GenomicSelection("hg19")
# chr1 and 2 from hg19, with a score column
GenomicSelection("hg19", c("chr1", "chr2"), "score")
}
\keyword{ manip }
