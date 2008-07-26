\name{hsTargets}
\alias{hsTargets}
\docType{data}
\title{ Human microRNAs and their target IDs}
\description{
 A set of human microRNA names and their corresponding known targets
 given as ensembl Trascript IDs.
}
\usage{data(hsTargets)}
\format{
 A data frame of microRNAs and their target ensembl IDs as recovered
 from miRBase.  Additional columns are also provided to  give the
 Chromosome as well as the start and end position of the microRNA
 binding site, and the strand orientation (plus or minus).
}
\details{
 Each mapping represents a different human microRNA, paired with one
 viable target.  Other information about where the microRNA binds is
 also included.  Some microRNAs have multiple targets and so some
 microRNAs may be represented more than once.
}
\source{
\url{http://microrna.sanger.ac.uk/sequences/index.shtml}
}
\references{
miRBase: microRNA sequences, targets and gene nomenclature.
Griffiths-Jones S, Grocock RJ, van Dongen S, Bateman A, Enright AJ.
NAR, 2006, 34, Database Issue, D140-D144

The microRNA Registry.
Griffiths-Jones S.
NAR, 2004, 32, Database Issue, D109-D111

}
\examples{
data(hsTargets)
}
\keyword{datasets}
