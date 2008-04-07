\name{laySequence-methods}
\docType{methods}
\alias{laySequence}
\alias{laySequence-methods}
\alias{laySequence,argoSession,character-method}
\title{ Load a sequence }
\description{
  Methods for loading sequences.
}
\section{Methods}{
\describe{
  The following methods are defined by \pkg{rtracklayer} for the
  \code{laySequence(object, sequence, ...)} generic.
  
  \item{object = "argoSession", sequence = "character"}{
    \code{laySequence(object, sequence, name, label = name)}: Load
    \code{sequence} into \code{object} under the name \code{name} and
    labeled by \code{label} in the user interface.
  }
}}
\keyword{methods}
