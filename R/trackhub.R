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

stopIfNotLocal <- function(x) {
    if (!uriIsWritable(x)) {
        stop("Repository is read only; cannot write on remote repository")
    }
}

hubFile <- function(x) paste(trimSlash(uri(x)), "hub.txt", sep = "/")
combineURI <- function(x,y) paste(trimSlash(uri(x)), y, sep = "/")

getGenomesContentList <- function(x) {
    if (uriExists(hubFile(x))) {
        hubContent <- getHubContent(hubFile(x))
        genomesFileValue <- hubContent["genomesFile"]
        if (!isFieldEmpty(genomesFileValue)) {
            genomesFilePath <- combineURI(x, genomesFileValue)
            genomesList <- getGenomesContent(genomesFilePath)
        }
        else message("hub.txt: 'genomesFile' does not contain valid reference to genomes file")
    }
}

setGenomesKey <- function(x, key, value) {
    genomesList <- getGenomesContentList(trackhub(x))
    position <- which(sapply(genomesList, function(y) genome(x) %in% y))
    genomesList[[position]][[key]] <- value
    hubContent <- getHubContent(hubFile(trackhub(x)))
    genomesFilePath <- combineURI(trackhub(x), hubContent["genomesFile"])
    cat("", file = genomesFilePath)
    sapply(genomesList, function(x) {
        setGenomeContentList(x, genomesFilePath)
    })
    return(invisible(NULL))
}

getGenomesKey <- function(x, key) {
    genomesList <- getGenomesContentList(trackhub(x))
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

setReplaceMethod("hub", "TrackHub", function(x, value) {
    stopIfNotLocal(hubFile(x))
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
    stopIfNotLocal(hubFile(x))
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
    stopIfNotLocal(hubFile(x))
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
    stopIfNotLocal(hubFile(x))
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
    stopIfNotLocal(hubFile(x))
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
    stopIfNotLocal(hubFile(x))
    hubContent <- getHubContent(hubFile(x))
    hubContent["descriptionUrl"] <- value
    setHubContent(x, hubContent)
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

setGeneric("addGenomesField", function(x, key, value) standardGeneric("addGenomesField"))
setGeneric("cols", function(x) standardGeneric("cols"))

setClass("TrackHubGenome",
         representation(trackhub = "TrackHub",
                        genome = "character"),
         contains = "TrackDb")

trackhub <- function(x) x@trackhub

getTrackDbContent <- function(x) {
    Tracks <- TrackContainer()
    contentDf <- readAndSanitize(x)
    contentDf$V1 <- gsub("^(\\t)+", "", contentDf$V1)
    tracksIndex <- grep("track", contentDf$V1)
    totalTracks <- length(tracksIndex)
    tracksIndex[[length(tracksIndex) + 1]] <- length(contentDf$V1) + 1 # to read last track from file
    trackNo <- 1L
    position <- 1L
    # to speed up, reading track by track
    while(trackNo <= totalTracks) {
        startPosition <- tracksIndex[trackNo]
        endPosition <- tracksIndex[trackNo + 1] - 1
        # dynamically creating track object
        track <- paste0(contentDf$V1[startPosition:endPosition],
                        "= '", contentDf$V2[startPosition:endPosition], "'", collapse=",")
        track <- paste0("Track(", track, ")")
        track <- eval(parse(text = track))
        Tracks[[position]] <- track
        position <- position + 1
        trackNo <- trackNo + 1
    }
    Tracks
}

createTrackHubGenome <- function(x, genomeRecord) {
    hubContent <- getHubContent(hubFile(trackhub(x)))
    genomesFilePath <- combineURI(trackhub(x), hubContent["genomesFile"])
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

setMethod("cols", "TrackHubGenome", function(x) {
    trackDbValue <- getGenomesKey(x, "trackDb")
    if (!isFieldEmpty(trackDbValue)) {
        trackDbFilePath <- combineURI(trackhub(x), trackDbValue)
        trackDbContent <- getTrackDbContent(trackDbFilePath)
    }
    else stop("genomes.txt: 'trackDb' does not contain valid reference to trackDb file")
})

setMethod("names", "TrackHubGenome", function(x) {
    trackDbContent <- cols(x)
    names(trackDbContent)
})

setMethod("trackNames", "TrackHubGenome", function(object) {
    names(object)
})

setMethod("addGenomesField", "TrackHubGenome", function(x, key, value) {
    genomesFields <- c("twoBitPath", "groups", "htmlPath", "metaDb", "trackDb" , "metaTab")
    if (key %in% genomesFields && !isFieldEmpty(value)) {
        createResource(combineURI(trackhub(x), value))
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
        twoBitFilePath <- combineURI(trackhub(x), twoBitPathValue)
        import(twoBitFilePath)
    }
    else stop("genome.txt: 'twoBitPath' does not contain a reference to a file")
})

setReplaceMethod("referenceSequence", "TrackHubGenome", function(x, value) {
    hubContent <- getHubContent(hubFile(trackhub(x)))
    genomesFilePath <- combineURI(trackhub(x), hubContent["genomesFile"])
    stopIfNotLocal(genomesFilePath)
    twoBitPathValue <- getGenomesKey(x, "twoBitPath")
    if (!isFieldEmpty(twoBitPathValue)) {
        twoBitFilePath <- combineURI(trackhub(x), twoBitPathValue)
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

setMethod("track", "TrackHubGenome", function(object, name, ...) {
    tracks <- cols(object)
    track <- sapply(tracks, function(x) {
        if(x@track == name)
            x
    })
    track <- Filter(Negate(is.null), track)
    if (length(track) == 1L) {
        if (isEmpty(track[[1]]@bigDataUrl)) {
            stop("Track '", name, "' does not contain any data file")
        }
        import(file.path(uri(object), track[[1]]@bigDataUrl))
    }else{
        stop("Track '", name, "' does not exist")
    }
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
                        priority = "character",
                        altColor = "character",
                        boxedCfg = "character",
                        chromosomes = "character",
                        darkerLabels = "character",
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
                        skipEmptyFields = "character",
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
                        indelDoubleInsert = "character",
                        indelQueryInsert = "character",
                        indelPolyA = "character",
                        minAliQual = "character",
                        pairEndsByName = "character",
                        pairSearchRange = "character",
                        showNames = "character",
                        doWiggle = "character",
                        maxWindowToDraw = "character",

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
                        itemRgb = "character",
                        colorByStrand = "character",
                        denseCoverage = "character",
                        labelOnFeature = "character",
                        exonArrows = "character",
                        exonNumbers = "character",
                        scoreFilter = "character",
                        scoreFilterLimits = "character",
                        maxItems = "character",
                        minGrayLevel = "character",
                        noScoreFilter = "character",
                        spectrum = "character",
                        scoreMax = "character",
                        scoreMin = "character",
                        thickDrawItem = "character",
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
                        exonArrowsDense = "character",
                        itemImagePath = "character",
                        itemBigImagePath = "character",
                        mergeSpannedItems = "character",
                        linkIdInName = "character",
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
                        showDiffBasesAllScales = "character",

                        # bigWig settings
                        autoscale = "character",
                        autoScale = "character",
                        viewLimits = "character",
                        viewLimitsMax = "character",
                        alwaysZero = "character",
                        graphTypeDefault = "character",
                        maxWindowToQuery = "character",
                        negateValues = "character",
                        smoothingWindow = "character",
                        transformFunc = "character",
                        windowingFunction = "character",
                        yLineMark = "character",
                        yLineOnOff = "character",
                        gridDefault = "character",


                        # hic settings
                        showSnpWidth = "character",
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
                        compositeTrack = "character",
                        allButtonPair = "character",
                        centerLabelsDense = "character",
                        dragAndDrop = "character",

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
                        viewUi = "character",
                        configurable = "character",

                        # multiWig
                        container = "character",
                        aggregate = "character",
                        showSubtrackColorOnUi = "character",

                        # Miscellaneous Deprecated Settings
                        metadata = "character",
                        noInherit = "character",
                        useScore = "character"
                        )
)

Track <- function(...) {
    new("Track", ...)
}

setClassUnion("Track_OR_TrackContainer", c("Track", "TrackContainer"))

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
