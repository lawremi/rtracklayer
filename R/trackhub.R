### =========================================================================
### TrackHub support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHub class
###

setClass("TrackHub", representation(uri = "character"),
         contains = "List")

uri <- function(x, ...) x@uri

getHubContent <- function(x)    {
  content <- readLines(x, warn = FALSE)
  rexp <- "^(\\w+)\\s?(.*)$"
  contentVec <- c(sub(rexp,"\\2", content))
  names(contentVec) <- sub(rexp, "\\1", content)
  contentVec
}

getGenomesContent <- function(x)    {
    content <- readLines(x, warn = FALSE)
    content_df <- read.csv(text = sub(" ", ",", content), header = FALSE)
    recordsList <- split(setNames(as.list(content_df$V2), content_df$V1),
                  cumsum(content_df$V1 == "genome"))
    recordsList
}

hubFile <- function(x) file.path(uri(x), "hub.txt")
genomesFile <- function(x, y) file.path(uri(x), y)

setMethod("genome", "TrackHub", function(x)    {
  hubContent <- getHubContent(hubFile(x))
  path <- genomesFile(x, hubContent["genomesFile"])
  genomeRecordsList <- getGenomesContent(path)
  genomes <- sapply(genomeRecordsList, function(x) x["genome"])
  as.character(genomes)
})

setMethod("[[", "TrackHub", function (x, i, j, ...)    {
  if (!missing(j))
    warning("argument 'j' ignored")
})

setMethod("$", "TrackHub", function (x, name)    {

})

setMethod("names", "TrackHub", function(x) genome(x))

setMethod("length", "TrackHub", function(x) length(names(x)))

setMethod("show", "TrackHub", function(object)    {
  cat(class(object), "repository\nuri:", uri(object), "\n")
  cat(S4Vectors:::labeledLine("genomes", genome(object)))
})

TrackHub <- function(uri, create = FALSE)    {
  if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  if (create)    {
    if (uriExists(uri))    {
      message("NOTE: '", uri, "' already exists")
      create <- FALSE
    } ## must create this before calling normURI (requires existence)
    else createResource(uri, dir = TRUE)
  }
  th <- new("TrackHub", uri = normURI(uri))
  if (create)
    createResource(hubFile(th))
  th
}

setAs("character", "TrackHub", function(from) TrackHub(from))
