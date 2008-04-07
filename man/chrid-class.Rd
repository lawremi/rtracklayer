\name{chrid-class}
\docType{class}
\alias{chrid-class}

\title{Class "chrid"}
\description{ An object representing a chromosome identifier. Meant to
  be used for subsetting \code{\linkS4class{trackSet}} instances by chromosome.}
\section{Objects from the Class}{
  Instances may be created through a call of the form \code{\link{chrid}(id)},
  where \code{id} is the string identifier of the chromosome (e.g. "chr22").
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} holding the
      string representation of the chromosome ID. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
Class \code{"\link[Biobase:characterORMIAME-class]{characterORMIAME}"}, by class "character", distance 2.
}
\section{Methods}{
  No methods defined with class "chrid" in the signature.
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{chrid}} for obtaining instances of this class.
  
  \code{\link{[.trackSet}} where instances of this class may be used for
  subsetting \code{\linkS4class{trackSet}} instances by chromosome.
}
\keyword{classes}
