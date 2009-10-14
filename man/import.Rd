\name{import}
\alias{import}
\alias{import,connection,character,ANY-method}
\alias{import,character,missing,ANY-method}
\alias{import,character,character,ANY-method}
\alias{import,missing,ANY,character-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Importing objects }
\description{
  Imports an object from a connection according to a specified format.
}
\usage{
import(con, format, text, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{ The connection through which the data is received. If this
    is a character vector, it is assumed to be a filename. }
  \item{format}{ The format in which to expect the input. If omitted and
    \code{con} is a filename, the format is taken from the file extension.}
  \item{text}{ If \code{con} is missing, this can be a character vector
    directly providing the string data to import. }
  \item{\dots}{ Arguments to pass to the format-specific import routines. }
}
\details{
  This function delegates to a format-specific function named according
  to the scheme \code{import.format} where \code{format} is specified by the
  \code{format} parameter.
}
\value{
  The object parsed from the connection or text.
}
\author{ Michael Lawrence }
\seealso{ \code{\link{export}} to do the reverse. }
\examples{
  track <- import(system.file("tests", "bed.wig", package = "rtracklayer"))
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"), version = "1")
  # or
  track <- import(system.file("tests", "v1.gff", package = "rtracklayer"), "gff1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
