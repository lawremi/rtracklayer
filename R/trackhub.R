### =========================================================================
### TrackHub support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHub class
###

setGeneric("uri", function(x) standardGeneric("uri"))

setClass("TrackHub", representation(uri = "character"),
         contains = "List")

setMethod("uri", "TrackHub", function(x) {
    x@uri
})


getHubContent <- function(x) {
    content <- readLines(x, warn = FALSE)
    rexp <- "^(\\w+)\\s?(.*)$"
    contentVec <- c(sub(rexp, "\\2", content))
    names(contentVec) <- sub(rexp, "\\1", content)
    contentVec
}

getGenomesContent <- function(x) {
    content <- readLines(x, warn = FALSE)
    content_df <- read.csv(text = sub(" ", ",", content), header = FALSE)
    recordsList <- split(setNames(as.list(content_df$V2), content_df$V1),
                         cumsum(content_df$V1 == "genome"))
    recordsList
}

isFileReference <- function(x) {
    tools::file_ext(x) == "txt"
}

hubFile <- function(x) file.path(uri(x), "hub.txt")
genomesFile <- function(x, y) file.path(uri(x), y)

setMethod("genome", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    genomesFileValue <- hubContent["genomesFile"]
    if ((isFileReference(genomesFileValue) && !is.na(genomesFileValue))) {
        genomesFilePath <- genomesFile(x, genomesFileValue)
        genomesRecordsList <- getGenomesContent(genomesFilePath)
        genomes <- sapply(genomesRecordsList, function(x) x["genome"])
        as.character(genomes)
    }
    else stop("hub.txt: 'genomesFile' does not contain valid reference to genomes file")
})

setMethod("[[", "TrackHub", function (x, i, j, ...) {
    if (!missing(j))
        warning("argument 'j' ignored")
    TrackHubGenome(x, i, ...)
})

setMethod("$", "TrackHub", function (x, name) {
    TrackHubGenome(x, name)
})

setMethod("names", "TrackHub", function(x) genome(x))

setMethod("length", "TrackHub", function(x) length(names(x)))

setMethod("show", "TrackHub", function(object) {
    cat(class(object), "repository\nuri:", uri(object), "\n")
    cat(S4Vectors:::labeledLine("genomes", genome(object)))
})

TrackHub <- function(uri, create = FALSE) {
    if (!isTRUEorFALSE(create))
        stop("'create' must be TRUE or FALSE")
    if (create) {
        if (uriExists(uri)) {
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHubGenome class
###

setClass("TrackHubGenome",
         representation(trackhub = "TrackHub",
                        genome = "character"),
         contains = "TrackDb")

trackhub <- function(x, ...) x@trackhub
trackDbFile <- function(x,y) file.path(uri(x), y)


setMethod("genome", "TrackHubGenome", function(x) x@genome)

setMethod("uri", "TrackHubGenome", function(x)
          file.path(uri(trackhub(x)), genome(x)))

setMethod("names", "TrackHubGenome", function(x) {
    #TODO
    NULL
})

setMethod("length", "TrackHubGenome", function(x) {
    length(names(x))
})

setMethod("show", "TrackHubGenome", function(object) {
    cat(class(object), "track database\ngenome:", genome(object), "\ntrackhub:",
        uri(trackhub(object)), "\n")
    cat(S4Vectors:::labeledLine("names", names(object)))
})

TrackHubGenome <- function(trackhub, genome, create = FALSE) {
    if (!isTRUEorFALSE(create))
        stop("'create' must be TRUE or FALSE")
    trackhub <- as(trackhub, "TrackHub")
    thg <- new("TrackHubGenome", trackhub = trackhub, genome = genome)
    if (create) {
        #TODO
    }
    thg
}
