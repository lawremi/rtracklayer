\name{trackSet-methods}
\docType{methods}
\alias{trackSet}
\alias{trackSet-methods}
\alias{trackSet,ucscSession-method}
\alias{trackSet,argoSession-method}
\alias{trackSet,AnnotatedDataFrame-method}
\alias{trackSet,character-method}
\alias{trackSet,data.frame-method}
\alias{trackSet,IRanges-method}

\title{Get a trackSet}
\description{
  Methods for getting a \code{\linkS4class{trackSet}}.
}
\usage{trackSet(object, ...)}
\arguments{
  \item{object}{An object from which to obtain a \code{trackSet}.}
  \item{\dots}{Arguments passed on to methods.}
}
\section{Methods}{
  \describe{
    The following methods are defined by \pkg{rtracklayer}.

    \item{object = "AnnotatedDataFrame"}{
      \code{trackSet(object, dataVals = NULL, ...)}:
      A constructor that sets \code{object} directly as the
      \code{featureData} slot and \code{dataVals} into the
      \code{dataVals} slot, if non-\code{NULL}. The arguments in
      \code{...} should correspond to slots in the \code{trackSet} class.
    }
    \item{object = "data.frame"}{
      \code{trackSet(object, dataVals = NULL, ...)}: coerces
      \code{object} to \code{AnnotatedDataFrame} and delegates to it.
    }

    \item{object = "IRanges"}{
      \code{trackSet(object, chrom, strand = NA, dataVals = NULL, ...)}:
      Constructor that takes the start and end of each feature from
      \code{object}. The chromsome, strand, data values and other slots
      are specified by \code{chrom}, \code{strand}, \code{dataVals} and
      \code{...}, respectively.
    }
    
    \item{object = "character"}{
      \code{trackSet(object, start = 1, end = start + span - 1, strand = NA,
        span = 1, dataVals = NULL, breaks = NULL, ...)}: a constructor
      that takes \code{object} to be the chromosome
      values. \code{start}, \code{end}, \code{strand} and
      \code{dataVals} specify the starts, ends, strands and data
      values, respectively, of the features. \code{span}, if given,
      should be a scalar integer indicating widths of the features
      If \code{end} is specified with \code{span}, it must be length one
      and then the region from \code{start} to \code{end} is broken up
      into features of width \code{span}. If \code{breaks} is specified,
      \code{span} must be omitted and \code{start} and \code{end} must
      be scalars. Then the features are formed by breaking the region
      from \code{start} to \code{end} at the positions specified in
      \code{breaks}. Way too complicated for its own good.
    }
    
    \item{object = "ucscSession"}{
      \code{trackSet(object, name, segment = genomeSegment(object), format =
        "bed", table = NULL)}:
      Get a track from the session by its \code{name} in format
      \code{format}. Some built-in tracks have multiple series, each
      identified by a \code{table} name.
    }
    
    \item{object = "argoSession"}{
      \code{trackSet(object, name, segment = genomeSegment(object))}:
      Get a track from the session by its \code{name}.
    }
  }
}
\seealso{
  \code{\link{layTrack}} for loading a \code{trackSet} into a session.
}
\examples{
\dontrun{
session <- browseGenome()
track <- import(system.file("tests", "v1.gff", package = "rtracklayer"))
session <- layTrack(session, track, "My Track", format = "gff1", view = TRUE)
track <- trackSet(session, name = "My Track")
}
}
\keyword{methods}
