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

isFieldEmpty <- function(x) {
    if ((isFileReference(x) && !is.na(x))) {
        return(FALSE)
    }
    return(TRUE)
}

hubFile <- function(x) paste(uri(x), "hub.txt", sep = "/")
genomesFile <- function(x, y) paste(uri(x), y, sep = "/")

genomesContentList <- function(x) {
    if (uriExists(hubFile(x))) {
        hubContent <- getHubContent(hubFile(x))
        genomesFileValue <- hubContent["genomesFile"]
        if (!isFieldEmpty(genomesFileValue)) {
            genomesFilePath <- genomesFile(x, genomesFileValue)
            genomesList <- getGenomesContent(genomesFilePath)
        }
        else message("hub.txt: 'genomesFile' does not contain valid reference to genomes file")
    }
}

setMethod("genome", "TrackHub", function(x) {
    genomesList <- genomesContentList(x)
    genomes <- sapply(genomesList, function(x) x["genome"])
    as.character(genomes)
})

setMethod("getListElement", "TrackHub", function(x, i, exact = TRUE) {
    TrackHubGenome(x, i)
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
trackDbFile <- function(x,y) paste(uri(x), y, sep = "/")

setMethod("genome", "TrackHubGenome", function(x) x@genome)

setMethod("uri", "TrackHubGenome", function(x)
          paste(uri(trackhub(x)), genome(x), sep = "/"))

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
