\name{ucscTable}
\alias{ucscTable}
\alias{ucscTable,UCSCSession-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get a database table from UCSC }
\description{
  A generic function for getting a database table from
  \acronym{UCSC}. This is mostly for internal use though some may find
  it useful.
}
\usage{
ucscTable(object, range, track, table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The e.g. \code{\linkS4class{UCSCSession}} }
  \item{range}{ The \code{\linkS4class{RangesList}} for which to
    retrieve the data. }
  \item{track}{ The name of the track to which the table belongs. }
  \item{table}{ The name of the table. }
}
\value{
  A \code{data.frame} of the table.
}
\references{ The UCSC Table Browser:
  \url{http://genome.ucsc.edu/cgi-bin/hgTables}.
}
\author{ Michael Lawrence }
\seealso{ \code{\link{track}} for getting tracks,
  \code{\link{genomeSequence}} for getting sequences. }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interface }
