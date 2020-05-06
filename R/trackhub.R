### =========================================================================
### TrackHub support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHub class
###

setClass("TrackHub", representation(uri = "character"))

isFileReference <- function(s) {
  length(grep(".txt", s)) == 1
}

uri <- function(x, ...) x@uri

getContent <- function(x) {
  content <- readLines(x, warn = FALSE)
  rexp <- "^(\\w+)\\s?(.*)$"
  data.frame(field = sub(rexp, "\\1", content),
             value = sub(rexp,"\\2", content))
}

getValue <- function(field, content) {
  position <- grep(field, content$field)
  content$value[position]
}

hubFile <- function(x) file.path(uri(x), "hub.txt")

genomesFile <- function(x, y) file.path(uri(x), y)

setMethod("genome", "TrackHub", function(x) {
  hubContent <- getContent(hubFile(x))
  field <- "genomesFile"
  genomesFileValue <- getValue(field, hubContent)
  if (isFileReference(genomesFileValue)) {
    path <- genomesFile(x, genomesFileValue)
    genomeContent <- getContent(path)
    genomes <- getValue("genome", genomeContent)
  }
  else genomes <- genomesFileValue
  as.character(structure(genomes, name = field))
})

setMethod("length", "TrackHub", function(x) length(genome(x)))

setMethod("show", "TrackHub", function(object) {
  cat(class(object), "repository\nuri:", uri(object), "\n")
  cat(S4Vectors:::labeledLine("genomes", genome(object)))
})

TrackHub <- function(uri="TrackHub", create = FALSE){
   if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  if (create) {
    if (uriExists(uri)) {
      message("NOTE: '", uri, "' already exists")
      create <- FALSE
    } ## must create this before calling normURI (requires existence)
    else createResource(uri, dir = TRUE)
  }
  ql <- new("TrackHub", uri = normURI(uri))
  if (create)
    createResource(hubFile(ql))
  ql
}

setAs("character", "TrackHub", function(from) TrackHub(from))
