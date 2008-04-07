\name{genomeSegment-methods}
\docType{methods}
\alias{genomeSegment}
\alias{genomeSegment<-}
\alias{genomeSegment<-,ucscView-method}
\alias{genomeSegment<-,argoView-method}
\alias{genomeSegment-methods}
\alias{genomeSegment,character-method}
\alias{genomeSegment,missing-method}
\alias{genomeSegment,trackSet-method}
\alias{genomeSegment,trackSets-method}
\alias{genomeSegment,browserSession-method}
\alias{genomeSegment,ucscSession-method}
\alias{genomeSegment,ucscView-method}
\alias{genomeSegment,argoView-method}
\alias{genomeSegment,ucscCart-method}
\title{ Accessing a segment of a genome }
\description{
  Methods for creating, getting or setting a \code{\linkS4class{genomeSegment}}.
}
\section{Methods}{
  The following methods are defined by \pkg{rtracklayer}.
  
  \describe{
    
    \item{object = "character"}{
      \code{genomeSegment(object, chrom = character(0), start =
        character(0), end = character(0), segment = genomeSegment())}:
      Create a \code{genomeSegment} instance with the same slot values
      as \code{segment}, except where overriden by one of the other
      parameters, which all correspond to slots in
      \code{genomeSegment}. Note that \code{object} corresponds to the
      \code{genome} slot.
    }

    \item{object = "missing"}{
      \code{genomeSegment(object, genome = character(0), chrom =
        character(0), start = character(0), end = character(0),
        segment = new("genomeSegment"))}:
      Similar to above, except \code{object} is omitted and
      \code{genome} is an explicit parameter.
    }
    
    \item{object = "trackSet"}{ Get the segment spanned by a track. }
    
    \item{object = "trackSets"}{ Get the union (including gaps) of spans
      for all tracks. }
    
    \item{object = "browserSession"}{ Get the current segment
      (i.e. segment displayed by active view or default segment). }
    
    \item{object = "ucscSession"}{ Get the last accessed segment. }
    
    \item{object = "ucscView" or "argoView"}{
      \code{genomeSegment(object)}: Get the segment displayed by the
      view.
      
      \code{genomeSegment(object) <- value}: Set the segment displayed
      by the the view.
    }
}}
\keyword{methods}
\examples{
  ## Create a genome segment for genome "hg18", chromosome 22
  genomeSegment("hg18", "chr22")
  ## Explicitly specify parameters
  segment <- genomeSegment(genome = "hg18", chrom = "chr22", start = 150000)
  ## Add an 'end' value
  genomeSegment(end = 200000, segment = segment)
}
