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
setGeneric("genomesFile", function(x) standardGeneric("genomesFile"))
setGeneric("genomesFile<-", function(x, value) standardGeneric("genomesFile<-"))
setGeneric("email", function(x) standardGeneric("email"))
setGeneric("email<-", function(x, value) standardGeneric("email<-"))
setGeneric("descriptionUrl", function(x) standardGeneric("descriptionUrl"))
setGeneric("descriptionUrl<-", function(x, value) standardGeneric("descriptionUrl<-"))
setGeneric("writeTrackHub", function(x) standardGeneric("writeTrackHub"))

setClass("TrackHub",
         representation(
                        uri = "character",
                        hub = "character",
                        shortLabel = "character",
                        longLabel = "character",
                        genomesFile = "character",
                        email = "character",
                        descriptionUrl = "character"),
         prototype(
                   hub = NA_character_,
                   shortLabel = NA_character_,
                   longLabel = NA_character_,
                   genomesFile = NA_character_,
                   email = NA_character_,
                   descriptionUrl = NA_character_),
         contains = "List")

hubFile <- function(x) paste(trimSlash(uri(x)), "hub.txt", sep = "/")

stopIfNotLocal <- function(x) {
    if (!uriIsWritable(x)) {
        stop("Repository is read only; cannot write on remote repository")
    }
}

getHubContent <- function(x) {
    content <- readLines(hubFile(x), warn = FALSE)
    rexp <- "^(\\w+)\\s?(.*)$"
    contentVec <- c(sub(rexp, "\\2", content))
    names(contentVec) <- sub(rexp, "\\1", content)
    x@hub <- contentVec["hub"]
    x@shortLabel <- contentVec["shortLabel"]
    x@longLabel <- contentVec["longLabel"]
    x@genomesFile <- contentVec["genomesFile"]
    x@email <- contentVec["email"]
    x@descriptionUrl <- contentVec["descriptionUrl"]
    x
}

setHubContent <- function(x) {
    file = hubFile(x)
    cat("", file = file)
    if (!is.na(x@hub))
        cat("hub ", x@hub, "\n", append = TRUE, sep = "", file = file)
    if (!is.na(x@shortLabel))
        cat("shortLabel ", x@shortLabel, "\n", append = TRUE, sep = "", file = file)
    if (!is.na(x@longLabel))
        cat("longLabel ", x@longLabel, "\n", append = TRUE, sep = "", file = file)
    if (!isFieldEmpty(x@genomesFile))
        cat("genomesFile ", x@genomesFile, "\n", append = TRUE, sep = "", file = file)
    if (!is.na(x@email))
        cat("email ", x@email, "\n", append = TRUE, sep = "", file = file)
    if (!is.na(x@descriptionUrl))
        cat("descriptionUrl ", x@descriptionUrl, "\n", append = TRUE, sep = "", file = file)
}

getGenomesContentList <- function(x) {
    if (uriExists(hubFile(x))) {
        genomesFileValue <- x@genomesFile
        if (!isFieldEmpty(genomesFileValue)) {
            genomesFilePath <- combineURI(uri(x), unname(genomesFileValue))
            content <- readLines(genomesFilePath, warn = FALSE)
            content_df <- read.csv(text = sub(" ", ",", content), header = FALSE)
            genomesList <- split(setNames(as.list(content_df$V2), content_df$V1),
                         cumsum(content_df$V1 == "genome"))
        }
        else message("hub.txt: 'genomesFile' does not contain valid reference to genomes file")
    }
}

setGenomeContentList <- function(x, file) {
    sapply(names(x), function(y) {
        if(x[[y]] != "")
            cat(y, " ", x[[y]], "\n", append = TRUE, sep = "", file = file)
    })
    cat("\n", append = TRUE, file = file)
}

Genome <- function(genome = "", trackDb = "", twoBitPath = "", groups = "",
                   description = "", organism = "", defaultPos = "",
                   orderKey = "", scientificName = "", htmlPath = "") {
        c(as.list(environment()))
}

setMethod("uri", "TrackHub", function(x) {
    x@uri
})

setMethod("hub", "TrackHub", function(x) {
    unname(x@hub)
})

setReplaceMethod("hub", "TrackHub", function(x, value) {
    x@hub <- value
    x
})

setMethod("shortLabel", "TrackHub", function(x) {
    unname(x@shortLabel)
})

setReplaceMethod("shortLabel", "TrackHub", function(x, value) {
    x@shortLabel <- value
    x
})

setMethod("longLabel", "TrackHub", function(x) {
    unname(x@longLabel)
})

setReplaceMethod("longLabel", "TrackHub", function(x, value) {
    x@longLabel <- value
    x
})

setMethod("genomesFile", "TrackHub", function(x) {
    unname(x@genomesFile)
})

setReplaceMethod("genomesFile", "TrackHub", function(x, value) {
    x@genomesFile <- value
    x
})

setMethod("email", "TrackHub", function(x) {
    unname(x@email)
})

setReplaceMethod("email", "TrackHub", function(x, value) {
    x@email <- value
    x
})

setMethod("descriptionUrl", "TrackHub", function(x) {
    unname(x@descriptionUrl)
})

setReplaceMethod("descriptionUrl", "TrackHub", function(x, value) {
    x@descriptionUrl <- value
    x
})

setMethod("genome", "TrackHub", function(x) {
    genomesList <- getGenomesContentList(x)
    genomes <- sapply(genomesList, function(x) x["genome"])
    as.character(genomes)
})

setMethod("getListElement", "TrackHub", function(x, i, exact = TRUE) {
    TrackHubGenome(x, i)
})

setMethod("names", "TrackHub", function(x) genome(x))

setMethod("length", "TrackHub", function(x) length(names(x)))

setMethod("writeTrackHub", "TrackHub", function(x) {
    stopIfNotLocal(hubFile(x))
    setHubContent(x)
})

setMethod("show", "TrackHub", function(object) {
    cat(class(object), "repository\nuri:", uri(object), "\n")
    cat(S4Vectors:::labeledLine("genomes", genome(object)))
    cat("hub:", hub(object), "\n")
    cat("shortLabel:", shortLabel(object), "\n")
    cat("longLabel:", longLabel(object), "\n")
    cat("genomesFile:", genomesFile(object), "\n")
    cat("email:", email(object), "\n")
    cat("descriptionUrl:", descriptionUrl(object), "\n")
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
    th <- new("TrackHub")
    th@uri <- normURI(uri)
    if (create) {
        createResource(hubFile(th))
    }
    else th <- getHubContent(th)
    th
}

setAs("character", "TrackHub", function(from) TrackHub(from))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHubGenome class
###

setGeneric("setGenomesField", function(x, key, value) standardGeneric("setGenomesField"))
setGeneric("getTracks", function(x) standardGeneric("getTracks"))

setClass("TrackHubGenome",
         representation(trackhub = "TrackHub",
                        genome = "character",
                        tracks = "TrackContainer",
                        levels = "integer"),
         contains = "TrackDb")

trackhub <- function(x) x@trackhub

setGenomesKey <- function(x, key, value) {
    trackhub <- trackhub(x)
    genomesList <- getGenomesContentList(trackhub)
    position <- which(sapply(genomesList, function(y) genome(x) %in% y))
    genomesList[[position]][[key]] <- value
    genomesFilePath <- combineURI(uri(trackhub(x)), trackhub@genomesFile)
    cat("", file = genomesFilePath)
    sapply(genomesList, function(x) {
        setGenomeContentList(x, genomesFilePath)
    })
    return(invisible(NULL))
}

getGenomesKey <- function(x, key) {
    genomesList <- getGenomesContentList(trackhub(x))
    position <- which(sapply(genomesList, function(y) genome(x) %in% y))
    if (isEmpty(position))
        stop("Genome : '",genome(x), "' not found")
    genomesList[[position]][[key]]
}

createTrack <- function(trackDf) {
    fieldToType <- list(
        track = "character", type = "character", shortLabel = "character", longLabel = "character",
        bigDataUrl = "character", html = "character", visibility = "character", meta = "character",
        color = "character", priority = "numeric", altColor = "character", boxedCfg = "logical",
        chromosomes = "character", darkerLabels = "logical", dataVersion = "character",
        directUrl = "character", iframeUrl = "character", iframeOptions = "character",
        mouseOverField = "character", otherDb = "character", pennantIcon = "character",
        tableBrowser = "character", url = "character", urlLabel = "character", urls = "character",
        skipEmptyFields = "logical", skipFields = "character", sepFields = "character",
        refUrl = "character", bigDataIndex = "character", bamColorMode = "character",
        bamGrayMode = "character", aliQualRange = "character", baseQualRange = "character",
        bamColorTag = "character", noColorTag = "character", bamSkipPrintQualScore = "character",
        indelDoubleInsert = "logical", indelQueryInsert = "logical", indelPolyA = "logical",
        minAliQual = "character", pairEndsByName = "character", pairSearchRange = "character",
        showNames = "logical", doWiggle = "logical", maxWindowToDraw = "integer",
        barChartBars = "character", barChartColor = "character", barChartLabel = "character",
        barChartMaxSize = "character", barChartSizeWindows = "character", barChartMetric = "character",
        barChartUnit = "character", barChartMatrixUrl = "character", barChartSampleUrl = "character",
        maxLimit = "character", labelFields = "character", defaultLabelFields = "character",
        itemRgb = "logical", colorByStrand = "character", denseCoverage = "integer",
        labelOnFeature = "logical", exonArrows = "logical", exonNumbers = "logical",
        scoreFilter = "character", scoreFilterLimits = "character", maxItems = "integer",
        minGrayLevel = "character", noScoreFilter = "logical", spectrum = "logical",
        scoreMax = "integer", scoreMin = "integer", thickDrawItem = "logical", searchIndex = "character",
        searchTrix = "character", labelSeparator = "character", bedNameLabel = "character",
        exonArrowsDense = "logical", itemImagePath = "character", itemBigImagePath = "character",
        mergeSpannedItems = "logical", linkIdInName = "logical", nextExonText = "character",
        prevExonText = "character", scoreLabel = "character", showTopScorers = "character",
        linkDataUrl = "character", interactDirectional = "character", interactUp = "character",
        interactMultiRegion = "character", maxHeightPixels = "character", speciesOrder = "character",
        frames = "character", summary = "character", baseColorUseCds = "character",
        baseColorUseSequence = "character", baseColorDefault = "character",
        showDiffBasesAllScales = "logical", autoscale = "character", autoScale = "character",
        viewLimits = "character", viewLimitsMax = "character", alwaysZero = "logical",
        graphTypeDefault = "character", maxWindowToQuery = "integer", negateValues = "logical",
        smoothingWindow = "character", transformFunc = "character", windowingFunction = "character",
        yLineMark = "character", yLineOnOff = "logical", gridDefault = "logical",
        showSnpWidth = "integer", otherSpecies = "character", minQual = "character", minFreq = "character",
        hapClusterEnabled = "character", hapClusterColorBy = "character", hapClusterTreeAngle = "character",
        hapClusterHeight = "character", applyMinQual = "character", superTrack = "character",
        parent = "character", compositeTrack = "logical", allButtonPair = "logical",
        centerLabelsDense = "logical", dragAndDrop = "character",
        hideEmptySubtracks = "logical", hideEmptySubtracksMultiBedUrl = "character",
        hideEmptySubtracksSourcesUrl = "character", hideEmptySubtracksLabel = "character",
        subGroup1 = "character", subGroup2 = "character", subGroup3 = "character", subGroup4 = "character",
        subGroup5 = "character", subGroup6 = "character", subGroup7 = "character", subGroup8 = "character",
        subGroup9 = "character", subGroups = "character", dimensions = "character",
        filterComposite = "character", dimensionAchecked = "character", dimensionBchecked = "character",
        sortOrder = "character", view = "character", viewUi = "logical", configurable = "logical",
        container = "character", aggregate = "character", showSubtrackColorOnUi = "logical",
        metadata = "character", noInherit = "logical", useScore = "integer")
    trackDf$value <- gsub("[Oo]n", "TRUE", trackDf$value)
    trackDf$value <- gsub("[Oo]ff", "FALSE", trackDf$value)
    args <- Map(as, trackDf$value, fieldToType[trackDf$field])
    names(args) <- trackDf$field
    track <- do.call(Track, args)
    track
}

getTabCountList <- function(contentdf) {
    matches <- gregexpr("^(\\t)+", contentdf)
    tabCountList <- vapply(matches, attr, integer(1L), "match.length")
    tabCountList
}

readAndSanitize <- function(filepath) {
    fileContent <- readLines(filepath, warn = FALSE)
    fileContent <- gsub("^(\\t)*#(.)*", "", fileContent) # to avoid reading commented tracks
    fileContent <- gsub(",", ";", fileContent)
    contentDf <- read.csv(text = sub(" ", ",", fileContent), header = FALSE)
    contentDf$V2 <- gsub(";", ",", contentDf$V2)
    nonEmptyContent <- vapply(contentDf$V2, function(x) x!="", logical(1))
    contentDf <- contentDf[nonEmptyContent,]
    contentDf
}

getTrackDbContent <- function(x, trackDbFilePath) {
    Tracks <- TrackContainer()
    contentDf <- readAndSanitize(trackDbFilePath)
    tracksIndex <- grep("\\btrack\\b", contentDf$V1)
    levels <- getTabCountList(contentDf$V1)
    levels <- levels[tracksIndex]
    levels <- as.integer(gsub(-1, 0, levels))
    totalTracks <- length(tracksIndex)
    tracksIndex[length(tracksIndex) + 1] <- length(contentDf$V1) + 1 # to read last track from file
    contentDf$V1 <- gsub("^(\\t)+", "", contentDf$V1)
    trackNo <- 1L
    position <- 1L
    # to speed up, reading track by track
    while(trackNo <= totalTracks) {
        startPosition <- tracksIndex[trackNo]
        endPosition <- tracksIndex[trackNo + 1] - 1
        trackDf <- setNames(data.frame(contentDf$V1[startPosition:endPosition],
                                     contentDf$V2[startPosition:endPosition]),
                          c("field", "value"))
        track <- createTrack(trackDf)
        Tracks[[position]] <- track
        position <- position + 1
        trackNo <- trackNo + 1
    }
    x@tracks <- Tracks
    x@levels <- levels
    x
}

createTrackHubGenome <- function(x, genomeRecord) {
    trackhub <- trackhub(x)
    genomesFilePath <- combineURI(uri(trackhub), trackhub@genomesFile)
    if (uriExists(genomesFilePath) && genome(x) %in% genome(trackhub)) {
        message("NOTE: Genome '", genome(x), "' already exists")
        return ()
    }
    createResource(uri(x), dir = TRUE)
    if (!uriExists(genomesFilePath))
        createResource(genomesFilePath)
    setGenomeContentList(genomeRecord, genomesFilePath)
}

setMethod("genome", "TrackHubGenome", function(x) x@genome)

setMethod("uri", "TrackHubGenome", function(x)
          paste(trimSlash(uri(trackhub(x))), genome(x), sep = "/"))

setMethod("getTracks", "TrackHubGenome", function(x) {
    x@tracks
})

setMethod("names", "TrackHubGenome", function(x) {
    as.character(names(getTracks(x)))
})

setMethod("trackNames", "TrackHubGenome", function(object) {
    names(object)
})

setMethod("setGenomesField", "TrackHubGenome", function(x, key, value) {
    genomesFields <- c("twoBitPath", "groups", "htmlPath", "metaDb", "trackDb" , "metaTab")
    if (key %in% genomesFields && !isFieldEmpty(value)) {
        createResource(combineURI(uri(trackhub(x)), value))
    }
    setGenomesKey(x, key, value)
})

setMethod("organism", "TrackHubGenome", function(object) {
    organism <- getGenomesKey(object, "organism")
    as.character(organism)
})

setMethod("referenceSequence", "TrackHubGenome", function(x) {
    twoBitPathValue <- getGenomesKey(x, "twoBitPath")
    if (!isFieldEmpty(twoBitPathValue)) {
        twoBitFilePath <- combineURI(uri(trackhub(x)), twoBitPathValue)
        import(twoBitFilePath)
    }
    else stop("genome.txt: 'twoBitPath' does not contain a reference to a file")
})

setReplaceMethod("referenceSequence", "TrackHubGenome", function(x, value) {
    trackhub <- trackhub(x)
    genomesFilePath <- combineURI(uri(trackhub), trackhub@genomesFile)
    stopIfNotLocal(genomesFilePath)
    twoBitPathValue <- getGenomesKey(x, "twoBitPath")
    if (!isFieldEmpty(twoBitPathValue)) {
        twoBitFilePath <- combineURI(uri(trackhub(x)), twoBitPathValue)
        export.2bit(value, twoBitFilePath)
        x
    }
    else stop("genome.txt: 'twoBitPath' does not contain a reference to a file")
})

setMethod("length", "TrackHubGenome", function(x) {
    length(names(x))
})

setMethod("writeTrackHub", "TrackHubGenome", function(x) {
    trackDbValue <- getGenomesKey(x, "trackDb")
    trackDbFilePath <- combineURI(uri(trackhub(x)), trackDbValue)
    stopIfNotLocal(trackDbFilePath)
    tabStrings <- sapply(x@levels, function(y) {
        paste(replicate(y, "\t"), collapse = "")
    })
    slots <- slotNames(Track())
    i <- 1L
    tracks <- sapply(x@tracks, function(y) {
        track <- sapply(slots, function(slotName) {
            slotValue <- slot(y, slotName)
            if (!isEmpty(slotValue)) {
                paste0(tabStrings[i], slotName, " ", slotValue)
            }
            else NULL
        })
        track[length(track) + 1] <- ""
        i <<- i + 1
        track
    })
    tracks <- as.character(Filter(Negate(is.null), tracks))
    tracks <- gsub("TRUE", "on", tracks)
    tracks <- gsub("FALSE", "off", tracks)
    writeLines(tracks, trackDbFilePath)
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
    thg <- new("TrackHubGenome")
    thg@trackhub <- trackhub
    thg@genome <- genome
    if (create) {
        createTrackHubGenome(thg, genomeRecord)
    }
    trackDbValue <- getGenomesKey(thg, "trackDb")
    if (!isFieldEmpty(trackDbValue)) {
        trackDbFilePath <- combineURI(uri(trackhub(thg)), trackDbValue)
        if (!uriExists(trackDbFilePath)) {
            createResource(trackDbFilePath)
        }else if (file.size(parseURI(trackDbFilePath)$path) != 1) {
            thg <- getTrackDbContent(thg, trackDbFilePath)
        }
    }
    thg
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import of tracks from Track Hub
###

setMethod("track", "TrackHubGenome", function(object, name, ...) {
    track <- sapply(object@tracks, function(x) {
        if(x@track == name)
            x
    })
    track <- Filter(Negate(is.null), track)
    if (length(track) == 1L) {
        if (isEmpty(track[[1]]@bigDataUrl)) {
            stop("Track '", name, "' does not contain any data file")
        }
        if (uriIsWritable(track[[1]]@bigDataUrl)) {
            import(paste0(parseURI(uri(trackhub(object)))$path, "/", track[[1]]@bigDataUrl))
        }else {
            import(track[[1]]@bigDataUrl)
        }
    }else {
        stop("Track '", name, "' does not exist")
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export of tracks to Track Hub
###

copyResourceToTrackHub <- function(object, uri) {
    parsed_uri <- .parseURI(uri)
    if (parsed_uri$scheme == "")
        uri <- paste0("file://", uri)
    filename <- basename(uri)
    object_uri <- .parseURI(uri(object))
    if (uriIsLocal(object_uri)) {
        dest_file <- file.path(object_uri$path, filename)
        if (paste(uri(object), filename, sep = "/") != uri)
            ### FIXME: URLdecode() here because of R bug
            download.file(URLdecode(uri), dest_file)
    }
    else stop("TrackHub is not local; cannot copy track")
    filename
}

.exportToTrackHub <- function(object, name,
                              format = bestFileFormat(value, object),
                              index = TRUE, ..., value)
{
    filename <- paste(name, format, sep = ".")
    path <- paste(uri(object), filename, sep = "/")
    file <- export(value, path, format = format, index = index, ...)
    track(object, name, index = FALSE) <- file
    object
}

setReplaceMethod("track",signature(object = "TrackHubGenome", value = "ANY"),
                 .exportToTrackHub)

setReplaceMethod("track",
                 signature(object = "TrackHubGenome", value = "RsamtoolsFile"),
                 function(object, name, ..., value)
                 {
                     if (missing(name))
                         name <- basename(path(value))
                     track(object, name) <- URLencode(path(value))
                     copyResourceToTrackHub(object, URLencode(index(value)))
                     object
                 })

setReplaceMethod("track",
                 signature(object = "TrackHubGenome", value = "RTLFile"),
                 function(object, name, ..., value)
                 {
                     if (missing(name))
                         name <- basename(path(value))
                     track(object, name) <- URLencode(path(value))
                     object
                 })

setReplaceMethod("track",
                 signature(object = "TrackHubGenome", value = "character"),
                 function(object, name = basename(object), ..., value)
                 {
                     filename <- copyResourceToTrackHub(object, value)
                     track <- sapply(object@tracks, function(x) {
                                if (x@track == name)
                                    return(TRUE)
                                return(FALSE)
                     })
                     trackDbValue <- getGenomesKey(object, "trackDb")
                     trackDbFilePath <- combineURI(uri(trackhub(object)), trackDbValue)
                     bigDataUrlValue <-  sub(uri(trackhub(object)), "", uri(object))
                     bigDataUrlValue <- sub("^/", "", bigDataUrlValue)
                     bigDataUrlValue <- paste(bigDataUrlValue, filename, sep = "/")
                     trackPosition <- which(track)
                     trackDf <- setNames(data.frame(c("track", "bigDataUrl"),
                                                    c(name, bigDataUrlValue)),
                                         c("field", "value"))
                     track <- createTrack(trackDf)
                     if (isEmpty(trackPosition)) {
                         trackPosition <- length(object@tracks) + 1
                         object@tracks[[trackPosition]] <- track
                     }else {
                         object@tracks[[trackPosition]] <- track
                     }
                     object
                 })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackContainer class
###

setClass("TrackContainer",
         representation("SimpleList"),
         prototype(elementType = "Track_OR_TrackContainer")
)

setMethod("names", "TrackContainer", function(x) {
    vapply(x, function(y) y@track, "character" ,USE.NAMES = FALSE)
})

TrackContainer <- function() {
    new("TrackContainer")
}

setClass("Track",
         representation(
                        # common trackDb settings
                        track = "character",
                        type = "character",
                        shortLabel = "character",
                        longLabel = "character",
                        bigDataUrl = "character",
                        html = "character",
                        visibility = "character",
                        meta = "character",


                        # common optional settings
                        color = "character",
                        priority = "numeric",
                        altColor = "character",
                        boxedCfg = "logical",
                        chromosomes = "character",
                        darkerLabels = "logical",
                        dataVersion = "character",
                        directUrl = "character",
                        iframeUrl = "character",
                        iframeOptions = "character",
                        mouseOverField = "character",
                        otherDb = "character",
                        pennantIcon = "character",
                        tableBrowser = "character",
                        url = "character",
                        urlLabel = "character",
                        urls = "character",
                        skipEmptyFields = "logical",
                        skipFields = "character",
                        sepFields = "character",


                        ##settings by track type

                        # bam settings
                        refUrl = "character",
                        bigDataIndex = "character",
                        bamColorMode = "character",
                        bamGrayMode = "character",
                        aliQualRange = "character",
                        baseQualRange = "character",
                        bamColorTag = "character",
                        noColorTag = "character",
                        bamSkipPrintQualScore = "character",
                        indelDoubleInsert = "logical",
                        indelQueryInsert = "logical",
                        indelPolyA = "logical",
                        minAliQual = "character",
                        pairEndsByName = "character",
                        pairSearchRange = "character",
                        showNames = "logical",
                        doWiggle = "logical",
                        maxWindowToDraw = "integer",

                        # bigBarChart settings
                        barChartBars = "character",
                        barChartColor = "character",
                        barChartLabel = "character",
                        barChartMaxSize = "character",
                        barChartSizeWindows = "character",
                        barChartMetric = "character",
                        barChartUnit = "character",
                        barChartMatrixUrl = "character",
                        barChartSampleUrl = "character",
                        maxLimit = "character",
                        labelFields = "character",
                        defaultLabelFields = "character",
                        itemRgb = "logical",
                        colorByStrand = "character",
                        denseCoverage = "integer",
                        labelOnFeature = "logical",
                        exonArrows = "logical",
                        exonNumbers = "logical",
                        scoreFilter = "character",
                        scoreFilterLimits = "character",
                        maxItems = "integer",
                        minGrayLevel = "character",
                        noScoreFilter = "logical",
                        spectrum = "logical",
                        scoreMax = "integer",
                        scoreMin = "integer",
                        thickDrawItem = "logical",
                        searchIndex = "character",
                        searchTrix = "character",
                        labelSeparator = "character",
                        # UNSUPPORTED fields
                        # filter.<fieldName>
                        # filterByRange.<fieldName>
                        # filterLimits.<fieldName>
                        # filterText.<fieldName>
                        # filterType.<fieldName>
                        # filterValues.<fieldName>
                        # filterValuesDefault.<fieldName>
                        # filterType.<fieldName>
                        # filterLabel.<fieldName>
                        bedNameLabel = "character",
                        exonArrowsDense = "logical",
                        itemImagePath = "character",
                        itemBigImagePath = "character",
                        mergeSpannedItems = "logical",
                        linkIdInName = "logical",
                        nextExonText = "character",
                        prevExonText = "character",
                        scoreLabel = "character",
                        showTopScorers = "character",

                        # bigChain settings
                        linkDataUrl = "character",

                        # bigInteract settings
                        interactDirectional = "character",
                        interactUp = "character",
                        interactMultiRegion = "character",
                        maxHeightPixels = "character",
                        speciesOrder = "character",
                        frames = "character",
                        summary = "character",

                        # bigNarrowPeak settings
                        # UNSUPPORTED fields
                        #scoreFilter
                        #pValueFilter
                        #qValueFilter
                        #signalFilter
                        #<column>FilterLimits
                        #<column>FilterByRange

                        # bigPsl settings
                        baseColorUseCds = "character",
                        baseColorUseSequence = "character",
                        baseColorDefault = "character",
                        showDiffBasesAllScales = "logical",

                        # bigWig settings
                        autoscale = "character",
                        autoScale = "character",
                        viewLimits = "character",
                        viewLimitsMax = "character",
                        alwaysZero = "logical",
                        graphTypeDefault = "character",
                        maxWindowToQuery = "integer",
                        negateValues = "logical",
                        smoothingWindow = "character",
                        transformFunc = "character",
                        windowingFunction = "character",
                        yLineMark = "character",
                        yLineOnOff = "logical",
                        gridDefault = "logical",


                        # hic settings
                        showSnpWidth = "integer",
                        otherSpecies = "character",


                        # vcfTabix settings
                        minQual = "character",
                        minFreq = "character",
                        hapClusterEnabled = "character",
                        hapClusterColorBy = "character",
                        hapClusterTreeAngle = "character",
                        hapClusterHeight = "character",
                        applyMinQual = "character",


                        ##Grouping tracks into sets and hierarchies

                        # Supertrack
                        superTrack = "character",
                        parent = "character",

                        # Composite Tracks
                        compositeTrack = "logical",
                        allButtonPair = "logical",
                        centerLabelsDense = "logical",
                        dragAndDrop = "character",
                        hideEmptySubtracks = "logical",
                        hideEmptySubtracksMultiBedUrl = "character",
                        hideEmptySubtracksSourcesUrl = "character",
                        hideEmptySubtracksLabel = "character",


                        # Subgroups
                        subGroup1 = "character",
                        subGroup2 = "character",
                        subGroup3 = "character",
                        subGroup4 = "character",
                        subGroup5 = "character",
                        subGroup6 = "character",
                        subGroup7 = "character",
                        subGroup8 = "character",
                        subGroup9 = "character",
                        subGroups = "character",
                        dimensions = "character",
                        filterComposite = "character",
                        dimensionAchecked = "character",
                        dimensionBchecked = "character",
                        sortOrder = "character",

                        # Subgroups Views
                        view = "character",
                        viewUi = "logical",
                        configurable = "logical",

                        # multiWig
                        container = "character",
                        aggregate = "character",
                        showSubtrackColorOnUi = "logical",

                        # Miscellaneous Deprecated Settings
                        metadata = "character",
                        noInherit = "logical",
                        useScore = "integer"
                        )
)

Track <- function(...) {
    new("Track", ...)
}

setClassUnion("Track_OR_TrackContainer", c("Track", "TrackContainer"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

combineURI <- function(x,y) paste(trimSlash(x), y, sep = "/")

isFileReference <- function(x) {
    formats <- c("txt", "2bit", "html")
    tools::file_ext(x) %in% formats
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
