# main track data structure - based on code by VJ Carey in encoDnaseI package

# two primary components:
# 1) an optional vector of data values (numeric or factor)
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
                   dataVals = matrix(),
                   genome = "hg18",
                   ...)
          {
            .Object@genome <- genome
            dataVals <- assayData[["dataVals"]]
            featureData <- as(featureData, "AnnotatedDataFrame")
            stopifnot(nrow(dataVals) == nrow(pData(featureData)))
            phenoData <- as(phenoData, "AnnotatedDataFrame")
            stopifnot(ncol(dataVals) == nrow(pData(phenoData)))
            
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
setGeneric("featStart", function(object)standardGeneric("featStart"))
setMethod("featStart", "trackSet", function(object) {
 pData(featureData(object))$featStart
})

setGeneric("featEnd", function(object)standardGeneric("featEnd"))
setMethod("featEnd", "trackSet", function(object) {
 pData(featureData(object))$featEnd
})

## chromosome names (presumably a factor)
setGeneric("featChrom", function(object)standardGeneric("featChrom"))
setMethod("featChrom", "trackSet", function(object) {
  pData(featureData(object))$featChrom
})

## sequence strand (+/-/NA)
setGeneric("featStrand", function(object)standardGeneric("featStrand"))
setMethod("featStrand", "trackSet", function(object) {
  pData(featureData(object))$featStrand
})


# chromosome IDs - for subsetting by chromosome
# a character vector (factors are not convenient for comparison/subsetting)
setClass("chrid", contains="character")

# for extracting/coercing chromosome IDs
setGeneric("chrid", function(object) standardGeneric("chrid"))

# get chromosome field as a 'chrid'
setMethod("chrid", "trackSet", function(object) chrid(featChrom(object)))

# label a character vector as representing chromosome IDs
setMethod("chrid", "character", function(object) {
  object <- sub("chr", "", object) # a common prefix
  new("chrid", object)
})
# fallback, attempt to coerce to character vector
setMethod("chrid", "ANY", function(object) chrid(as.character(object)))

# FIXME: support subsetting by a genomeSegment instance
setMethod("[", "trackSet", function(x, i, j, ..., drop=FALSE) {
  if (missing(i)) i <- seq_along(featureNames(x))
  if (missing(j)) j <- seq_along(sampleNames(x))
  if (is(i, "chrid")) i <- which(chrid(x) %in% i)
  x@assayData <- assayDataNew(dataVals = dataVals(x)[i,j,drop=FALSE])
  x@featureData <- featureData(x)[i,]
  x@phenoData <- phenoData(x)[j,]
  x
})

# convert trackSet to a data.frame, adding 'seqMid' column
setGeneric("trackData", function(object, ...) standardGeneric("trackData"))
setMethod("trackData", "trackSet", 
  function(object) {
    Y <- dataVals(object)
    df <- pData(featureData(object))
    Xs <- df[,c("featStart", "featEnd")]
    X <- apply(data.matrix(Xs), 1, mean)
    df <- cbind(df, featMid = X)
    cbind(df, dataVals = Y)
  })

setAs("trackSet", "data.frame", function(from) trackData(from))

# get the genome segment spanned by the track
setMethod("genomeSegment", "trackSet", function(object) {
  co <- featChrom(object)
  if (!(all(co == co[1])))
    stop("data on multiple chromosomes present; rangeLocs not meaningful")
  cs <- pData(featureData(object))$featStart
  ce <- pData(featureData(object))$featEnd
  genomeSegment(genome = object@genome, chrom = as.character(co[1]),
                start = min(as.numeric(cs)), end = max(as.numeric(ce)))
})

# use a trackSet to subset a BString
#setMethod("views", c("BString", "trackSet"))

# 'trackSets' is a list of trackSet instances
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
          function(object) do.call("c", lapply(object, genomeSegment)))
