# main track data structure - based on code by VJ Carey in encoDnaseI package

# two primary components:
# 1) an optional matrix of data values (numeric or factor)
# 2) feature metadata, like location

setClass("trackSet",
         representation(genome = "character"),
         prototype(
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), trackSet="1.0.0"))),
         "eSet")


setMethod("initialize", "trackSet",
          function(.Object,
                   assayData = assayDataNew(
                     dataVals = as.matrix(dataVals), ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   dataVals = matrix(nrow = ifelse(missing(featureData), 0,
                                       nrow(pData(featureData))), ncol = 0),
                   genome = "hg18",
                   ...)
          {
            if (!missing(assayData) && !missing(dataVals))
              stop("'assayData' and 'dataVals' cannot both be specified")
            .Object@genome <- genome            
            callNextMethod(.Object,
                           assayData = assayData,
                           phenoData = phenoData,
                           featureData = featureData,
                           experimentData = experimentData,
                           annotation = annotation)
          })

### convenience accessors

## the actual data values
setGeneric("dataVals", function(object)standardGeneric("dataVals"))
setMethod("dataVals", "trackSet", function(object) {
  get("dataVals", assayData(object))
})

## start positions
#setGeneric("start", function(object)standardGeneric("start"))
setMethod("start", "trackSet", function(x) {
  pData(featureData(x))$start
})

#setGeneric("end", function(object)standardGeneric("end"))
setMethod("end", "trackSet", function(x) {
  pData(featureData(x))$end
})

## chromosome names (presumably a factor)
setMethod("chrom", "trackSet", function(object) {
  pData(featureData(object))$chrom
})

## sequence strand (+/-/NA)
setGeneric("strand", function(object, ...) standardGeneric("strand"))
setMethod("strand", "trackSet", function(object) {
  pData(featureData(object))$strand
})

setMethod("genome", "trackSet", function(object) object@genome)

# chromosome IDs - for subsetting by chromosome
# a character vector (factors are not convenient for comparison/subsetting)
setClass("chrid", contains="character")

# for extracting/coercing chromosome IDs
setGeneric("chrid", function(object) standardGeneric("chrid"))

# get chromosome field as a 'chrid'
setMethod("chrid", "trackSet", function(object) chrid(chrom(object)))

# label a character vector as representing chromosome IDs
setMethod("chrid", "character", function(object) {
  object <- sub("chr", "", object) # a common prefix
  new("chrid", object)
})
# fallback, attempt to coerce to character vector
setMethod("chrid", "ANY", function(object) chrid(as.character(object)))

setMethod("[", "trackSet", function(x, i, j, ..., drop=FALSE) {
  if (missing(i)) i <- seq_along(featureNames(x))
  if (missing(j)) j <- seq_along(sampleNames(x))
  if (is(i, "chrid")) i <- which(chrid(x) %in% i)
  if (is(i, "genomeSegment")) {
    if (length(genome(i)) && genome(i) != genome(x))
      stop("Segment genome (", genome(i), ") and trackSet genome (", genome(x),
           ") do not match")
    i <- rep(TRUE, seq_along(featureNames(x)))
    if (length(chrom(i)))
      i <- i & (chrom(x) == chrom(i))
    if (length(start(i)))
      i <- i & (start(x) >= start(i))
    if (length(end(i)))
      i <- i & (end(x) <= end(i))
  }
  x@assayData <- assayDataNew(dataVals = dataVals(x)[i,j,drop=FALSE])
  x@featureData <- featureData(x)[i,]
  x@phenoData <- phenoData(x)[j,]
  x
})

setGeneric("split", function(x, f, drop = FALSE, ...) standardGeneric("split"))
setMethod("split", "trackSet", function(x, f, drop = FALSE) {
  splitInd <- split(seq_along(featureNames(x)), f)
  trackSets(lapply(splitInd, function(ind) x[ind,]))
})

# convert trackSet to a data.frame, adding 'featMid' column
setGeneric("trackData", function(object, ...) standardGeneric("trackData"))
setMethod("trackData", "trackSet", 
  function(object) {
    Y <- dataVals(object)
    df <- pData(featureData(object))
    Xs <- df[,c("start", "end")]
    X <- apply(data.matrix(Xs), 1, mean)
    df <- cbind(df, featMid = X)
    cbind(df, dataVals = Y)
  })

setAs("trackSet", "data.frame", function(from) trackData(from))

### convenience methods for getting and building AnnotatedDataFrames that
### conform to the requirements on the featureData of a trackSet
### not often used directly, see trackSet()

setGeneric("trackFeatureData",
           function(object, ...) standardGeneric("trackFeatureData"))

## extract the feature data with only the canonical location fields
setMethod("trackFeatureData", "trackSet",
          function(object)
          {
            fields <- c("chrom", "start", "end", "strand")
            featureData(object)[,fields]
          })

## annotate a data.frame for use as featureData in a trackSet
setMethod("trackFeatureData", "data.frame",
          function(object)
          {
            required <- c("chrom", "start", "end")
            if (!(all(required %in% colnames(object))))
              stop("Data frame must at least include columns: ",
                   paste(required, collapse = ", "), ".")
            if (is.null(object$strand))
              object$strand <- NA
            object$chrom <- paste("chr", sub("chr", "", object$chrom), sep = "")
            object$start <- as.integer(as.character(object$start))
            object$end <- as.integer(as.character(object$end))
            features <- as(object, "AnnotatedDataFrame")
            descs <- c(chrom = "Chromosome ID",
                       start = "Start position on chromosome",
                       end = "End position on chromosome",
                       strand = "DNA strand, sense (+) or antisense (-)")
            md <- data.frame(labelDescription = descs, stringsAsFactors = FALSE)
            varMetadata(features)[rownames(md),] <- md
            features
          })

## convenience for creating a featureData for a trackSet
## specify:
##  start and end for arbitrary regions (i.e. genes)
##  start and span for regions of given widths at different starts
##    (i.e. individual bases = 1, codons = 3)
##  end (length 1) and span (or one of them and dataVals), possibly start:
#     fixed-width splitting (i.e. a value for each base (= 1) or codon (= 3))
##  breaks and end, possibly start: segmentation points
setMethod("trackFeatureData", "character",
          function(object, start = 1, end = start + span - 1, strand = NA,
                   span = 1, dataVals = NULL, breaks = NULL, ...)
          {
            if (!missing(end) && !missing(span) || !is.null(dataVals)) {
              args <- list(from = start)
              if (!missing(end))
                args <- c(args, to = end)
              if (is.null(dataVals) || missing(end))
                args <- c(args, by = span)
              if (!is.null(dataVals))
                args <- c(args, along.with = dataVals)
              start <- do.call("seq", args)
              if (missing(span))
                end <- c(tail(start, -1) - 1, end)
              else end <- start + span - 1
            } else if (!is.null(breaks)) {
              start <- c(start, breaks)
              end <- c(tail(start, -1) - 1, end)
            }
            df <- data.frame(chrom = object, start = start,
                             end = end, strand = strand, ...)
            trackFeatureData(df)
          })

## or from an IRanges object
setMethod("trackFeatureData", "IRanges",
          function(object, chrom, strand = NA, ...)
          {
            trackFeatureData(chrom, start(object), end(object), strand,
                             ...)
          })

## get the genome segment spanned by the track
setMethod("genomeSegment", "trackSet", function(object) {
  segment <- genomeSegment(genome = genome(object))
  co <- chrom(object)
  if ((all(co == co[1]))) {
    cs <- start(object)
    ce <- end(object)
    segment <- genomeSegment(chrom = as.character(co[1]),
                             start = min(as.numeric(cs)),
                             end = max(as.numeric(ce)),
                             segment = segment)
  }
  segment
})

## get a trackSet
setGeneric("trackSet", function(object, ...) standardGeneric("trackSet"))

## conversion from an eSet
setMethod("trackSet", "eSet",
          function(object, span, dataVals = NULL, ...)
          {
            ## figure out chrom,start,end,strand for each record
            ann <- annotation(object)
            keys <- featureNames(object)
            locMap <- getAnnMap("CHRLOC", ann)[keys]
            locList <- as.list(locMap)
            locInd <- rep(seq_along(locList), sapply(locList, length))
            start <- unlist(locList)
            mapped <- !is.na(start)
            start <- start[mapped]
            chrom <- sub(".*?\\.([a-zA-Z0-9]*).*", "\\1", names(start))
            chrom <- paste("chr", chrom, sep = "")
            strand <- ifelse(start > 0, "+", "-")
            ## combine with existing featureData
            features <- trackFeatureData(chrom, abs(start), strand = strand,
                                         span = span)
            origFeatureData <- featureData(object)[locInd,][mapped,] 
            dimLabels(features) <- dimLabels(origFeatureData) # match dimLabels
            featureNames(features) <- featureNames(origFeatureData)
            featureData <- combine(origFeatureData, features)
            args <- list(featureData = featureData, genome = genome)
            ## extract data to use as 'dataVals'
            if (!is.null(dataVals)) {
              dataVals <- as.matrix(dataVals)[locInd,,drop=FALSE][mapped,]
              args <- c(args, list(dataVals = dataVals))
            }
            ## create new trackSet object
            do.call("new", c("trackSet", args))
          })

## ExpressionSet has different parameter defaults
setMethod("trackSet", "ExpressionSet",
          function(object, span = 25, dataVals = NULL, ...)
          {
            callNextMethod(object, span, dataVals, ...)
          })

## Some simple wrappers around trackFeatureData()
setMethod("trackSet", "data.frame",
          function(object, dataVals = NULL, ...) {
            trackSet(trackFeatureData(object), dataVals, ...)
          })
setMethod("trackSet", "character",
          function(object, start = 1, end = start + span - 1, strand = NA,
                   span = 1, dataVals = NULL, breaks = NULL, ...)
          {
            fd <- trackFeatureData(object, start, end, strand, span, dataVals,
                                   breaks)
            trackSet(fd, dataVals, ...)
          })
setMethod("trackSet", "IRanges",
          function(object, chrom, strand = NA, dataVals = NULL, ...)
          {
            fd <- trackFeatureData(chrom, start(object), end(object), strand)
            trackSet(fd, dataVals, ...)
          })

## core wrapper
setMethod("trackSet", "AnnotatedDataFrame",
          function(object, dataVals = NULL, ...)
          {
            args <- list(featureData = object, ...)
            if (!is.null(dataVals))
              args <- c(list(dataVals = dataVals), args)
            do.call("new", c("trackSet", args))
          })

## use a trackSet to create a BStringViews
setGeneric("trackViews",
           function(object, track, ...) standardGeneric("trackViews"))
setMethod("trackViews", c("XString", "trackSet"),
          function(object, track)
          {
            views(object, start(track), end(track))
          })

### show method
setMethod("show", "trackSet",
          function(object)
          {
            callNextMethod()
            cat("genome:", genome(object), "\n")
          })

### 'trackSets' is a list of trackSet instances
setClass("trackSets", contains = "list")

setGeneric("trackSets", function(object, ...) standardGeneric("trackSets"))

setMethod("trackSets", "trackSet",
          function(object, ...) new("trackSets", list(object, ...)))
setMethod("trackSets", "missing", function(object) trackSets(list()))
setMethod("trackSets", "list",
          function(object)
          {
            stopifnot(all(sapply(object, is, "trackSet")))
            new("trackSets", object)
          })

setMethod("genomeSegment", "trackSets",
          function(object) {
            if (!length(object))
              genomeSegment()
            else do.call("c", lapply(object, genomeSegment))
          })
