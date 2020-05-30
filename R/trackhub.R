### =========================================================================
### TrackHub support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHub class
###

setGeneric("hub", function(x) standardGeneric("hub"))
setGeneric("hub<-", function(x, value) standardGeneric("hub<-"))
setGeneric("shortLabel", function(x) standardGeneric("shortLabel"))
setGeneric("shortLabel<-", function(x, value) standardGeneric("shortLabel<-"))
setGeneric("longLabel", function(x) standardGeneric("longLabel"))
setGeneric("longLabel<-", function(x, value) standardGeneric("longLabel<-"))
setGeneric("genomeFile", function(x) standardGeneric("genomeFile"))
setGeneric("genomeFile<-", function(x, value) standardGeneric("genomeFile<-"))
setGeneric("email", function(x) standardGeneric("email"))
setGeneric("email<-", function(x, value) standardGeneric("email<-"))
setGeneric("descriptionUrl", function(x) standardGeneric("descriptionUrl"))
setGeneric("descriptionUrl<-", function(x, value) standardGeneric("descriptionUrl<-"))

setClass("TrackHub", representation(uri = "character"),
         contains = "List")

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
    tools::file_ext(x) == "txt" || tools::file_ext(x) == "2bit"
}

isFieldEmpty <- function(x) {
    if ((isFileReference(x) && !is.na(x)) && !is.null(x)) {
        return(FALSE)
    }
    return(TRUE)
}

trimSlash <- function(x) {
    sub("/$", "", x)
}

hubFile <- function(x) paste(trimSlash(uri(x)), "hub.txt", sep = "/")
genomesFile <- function(x, y) paste(trimSlash(uri(x)), y, sep = "/")

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

getGenomesKey <- function(x, key) {
    genomesList <- genomesContentList(trackhub(x))
    position <- which(sapply(genomesList, function(y) genome(x) %in% y))
    genomesList[[position]][[key]]
}

Genome <- function(genome = "", trackDb = "", twoBitPath = "", groups = "",
                   description = "", organism = "", defaultPos = "",
                   orderKey = "", scientificName = "", htmlPath = "") {
        c(as.list(environment()))
}

setGenomeContentList <- function(x, file) {
    sapply(names(x), function(y) {
        if(x[[y]] != "")
            cat(y, " ", x[[y]], "\n", append = TRUE, sep = "", file = file)
    })
    cat("\n", append = TRUE, file = file)
}

setHubContent <- function(x, hubContent) {
    cat("", file = hubFile(x))
    sapply(names(hubContent), function(y) {
        if (!is.null(hubContent[y]) && !is.na(hubContent[y]))
            cat(y, " ", hubContent[y], "\n", append = TRUE, sep = "", file = hubFile(x))
    })
    cat("\n", append = TRUE, file = hubFile(x))
}

setMethod("uri", "TrackHub", function(x) {
    x@uri
})

setMethod("hub", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    hubContent[["hub"]]
})

setReplaceMethod("hub", "TrackHub", function(x, value){
    hubContent <- getHubContent(hubFile(x))
    hubContent["hub"] <- value
    setHubContent(x, hubContent)
    x
})

setMethod("shortLabel", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    hubContent[["shortLabel"]]
})

setReplaceMethod("shortLabel", "TrackHub", function(x, value) {
    hubContent <- getHubContent(hubFile(x))
    hubContent["shortLabel"] <- value
    setHubContent(x, hubContent)
    x
})

setMethod("longLabel", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    hubContent[["longLabel"]]
})

setReplaceMethod("longLabel", "TrackHub", function(x, value) {
    hubContent <- getHubContent(hubFile(x))
    hubContent["longLabel"] <- value
    setHubContent(x, hubContent)
    x
})

setMethod("genomeFile", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    hubContent[["genomesFile"]]
})

setReplaceMethod("genomeFile", "TrackHub", function(x, value) {
    hubContent <- getHubContent(hubFile(x))
    hubContent["genomesFile"] <- value
    setHubContent(x, hubContent)
    x
})

setMethod("email", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    hubContent[["email"]]
})

setReplaceMethod("email", "TrackHub", function(x, value) {
    hubContent <- getHubContent(hubFile(x))
    hubContent["email"] <- value
    setHubContent(x, hubContent)
    x
})

setMethod("descriptionUrl", "TrackHub", function(x) {
    hubContent <- getHubContent(hubFile(x))
    hubContent[["descriptionUrl"]]
})

setReplaceMethod("descriptionUrl", "TrackHub", function(x, value) {
    hubContent <- getHubContent(hubFile(x))
    hubContent["descriptionUrl"] <- value
    setHubContent(x, hubContent)
    x
})

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
trackDbFile <- function(x,y) paste(trimSlash(uri(x)), y, sep = "/")
twobitFile <- function(x,y) paste(trimSlash(uri(x)), y, sep = "/")

getTrackDbContent <- function(x) {
    content <- readLines(x, warn = FALSE)
    # TODO
    # store content in right data structure
}

createTrackHubGenome <- function(x, genomeRecord) {
    hubContent <- getHubContent(hubFile(trackhub(x)))
    genomesFilePath <- genomesFile(trackhub(x), hubContent["genomesFile"])
    if (uriExists(genomesFilePath) && genome(x) %in% genome(trackhub(x))) {
        message("NOTE: Genome '", genome(x), "' already exists")
        return ()
    }
    createResource(uri(x), dir = TRUE)
    createResource(genomesFilePath)
    setGenomeContentList(genomeRecord, genomesFilePath)
}

setMethod("genome", "TrackHubGenome", function(x) x@genome)

setMethod("uri", "TrackHubGenome", function(x)
          paste(trimSlash(uri(trackhub(x))), genome(x), sep = "/"))

setMethod("names", "TrackHubGenome", function(x) {
    trackDbValue <- getGenomesKey(x, "trackDb")
    if (!isFieldEmpty(trackDbValue)) {
        trackDbFilePath <- trackDbFile(trackhub(x), trackDbValue)
        getTrackDbContent(trackDbFilePath)
    }
    else stop("genomes.txt: 'trackDb' does not contain valid reference to trackDb file")
})

setMethod("organism", "TrackHubGenome", function(object) {
    organism <- getGenomesKey(object, "organism")
    as.character(organism)
})

setMethod("referenceSequence", "TrackHubGenome", function(x) {
    twoBitPathValue <- getGenomesKey(x, "twoBitPath")
    if (!isFieldEmpty(twoBitPathValue)) {
        twoBitFilePath <- twobitFile(trackhub(x), twoBitPathValue)
        import(twoBitFilePath)
    }
    else stop("genome.txt: 'twoBitPath' does not contain a reference to a file")
})

setReplaceMethod("referenceSequence", "TrackHubGenome", function(x, value) {
    twoBitPathValue <- getGenomesKey(x, "twoBitPath")
    if (!isFieldEmpty(twoBitPathValue)) {
        twoBitFilePath <- twobitFile(trackhub(x), twoBitPathValue)
        export.2bit(value, twoBitFilePath)
        x
    }
    else stop("genome.txt: 'twoBitPath' does not contain a reference to a file")
})

setMethod("length", "TrackHubGenome", function(x) {
    length(names(x))
})

setMethod("show", "TrackHubGenome", function(object) {
    cat(class(object), "track database\ngenome:", genome(object), "\ntrackhub:",
        uri(trackhub(object)), "\n")
    cat(S4Vectors:::labeledLine("names", names(object)))
})

TrackHubGenome <- function(trackhub, genome, create = FALSE, genomeRecord = Genome()) {
    if (!isTRUEorFALSE(create))
        stop("'create' must be TRUE or FALSE")
    trackhub <- as(trackhub, "TrackHub")
    thg <- new("TrackHubGenome", trackhub = trackhub, genome = genome)
    if (create) {
        createTrackHubGenome(thg, genomeRecord)
    }
    thg
}
