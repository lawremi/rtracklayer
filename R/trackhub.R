### =========================================================================
### TrackHub support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackContainer class
###

setClass("TrackContainer",
         representation("SimpleList"),
         prototype(elementType = "Track")
)

setMethod("names", "TrackContainer", function(x) {
    vapply(x, function(y) y@track, character(1L) ,USE.NAMES = FALSE)
})

TrackContainer <- function(...) {
    args <- list(...)
    if (length(args) == 1 && is.list(args[[1L]]))
        args <- args[[1L]]
    if (!all(vapply(args, is, logical(1L), "Track")))
        stop("all elements in '...' must be Track objects")
    S4Vectors:::new_SimpleList_from_list("TrackContainer", args)
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Genome class
###

setClass("Genome",
         representation(
                        genome = "character",
                        trackDb = "character",
                        metaDb = "character",
                        metaTab = "character",
                        twoBitPath = "character",
                        groups = "character",
                        description = "character",
                        organism = "character",
                        defaultPos = "character",
                        orderKey = "character",
                        htmlPath = "character"
         ),
         prototype(
                   genome = NA_character_,
                   trackDb = NA_character_,
                   metaDb = NA_character_,
                   metaTab = NA_character_,
                   twoBitPath = NA_character_,
                   groups = NA_character_,
                   description = NA_character_,
                   organism = NA_character_,
                   defaultPos = NA_character_,
                   orderKey = NA_character_,
                   htmlPath = NA_character_
         )
)

Genome <- function(...) {
        new("Genome", ...)
}

stopIfNotGenome <- function(x) {
    if (!is(x, "Genome"))
        stop("value must be Genome object")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomeContainer class
###

setClass("GenomeContainer",
         representation("SimpleList"),
         prototype(elementType = "Genome")
)

setMethod("names", "GenomeContainer", function(x) {
    vapply(x, function(y) y@genome, character(1L), USE.NAMES = FALSE)
})

setMethod("getListElement", "GenomeContainer", function(x, i, exact = TRUE) {
    genome <- x[names(x) == i]
    if (length(genome) == 1L) unlist(genome)[[1L]]
})

GenomeContainer <- function(...) {
    args <- list(...)
    if (length(args) == 1 && is.list(args[[1L]]))
        args <- args[[1L]]
    if (!all(vapply(args, is, logical(1L), "Genome")))
        stop("all elements in '...' must be Genome objects")
    S4Vectors:::new_SimpleList_from_list("GenomeContainer", args)
}

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
setGeneric("genomeField", function(x, name, key) standardGeneric("genomeField"))
setGeneric("genomeField<-", function(x, name, key, value) standardGeneric("genomeField<-"))
setGeneric("genomeInfo", function(x, name) standardGeneric("genomeInfo"))
setGeneric("genomeInfo<-", function(x, value) standardGeneric("genomeInfo<-"))

setClass("TrackHub",
         representation(
                        uri = "character",
                        hub = "character",
                        shortLabel = "character",
                        longLabel = "character",
                        genomesFile = "character",
                        email = "character",
                        descriptionUrl = "character",
                        genomeContainer = "GenomeContainer"),
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

getGenomesContent <- function(x) {
    if (uriExists(hubFile(x))) {
        genomesFileValue <- x@genomesFile
        if (!isFieldEmpty(genomesFileValue)) {
            genomesFilePath <- combineURI(uri(x), unname(genomesFileValue))
            if (isFileEmpty(genomesFilePath)) {
                return(list())
            }
            content <- readLines(genomesFilePath, warn = FALSE)
            content_df <- read.csv(text = sub(" ", ",", content), header = FALSE)
            genomesIndex <- grep("\\bgenome\\b", content_df$V1)
            totalGenomes <- length(genomesIndex)
            genomesIndex[length(genomesIndex) + 1] <- length(content_df$V1) + 1
            genomes <- lapply(1:totalGenomes, function(x) {
                startPosition <- genomesIndex[x]
                endPosition <- genomesIndex[x + 1] - 1
                genome <- setNames(data.frame(content_df$V1[startPosition:endPosition],
                                    content_df$V2[startPosition:endPosition]),
                         c("field", "value"))
                genome <- setNames(as.list(genome$value), genome$field)
                genome <- do.call(Genome, genome)
            })
            genomes
        }
        else message("hub.txt: 'genomesFile' does not contain valid reference to genomes file")
    }
}

setGenomesContent <- function(x, genomeContainer) {
    genomesFilePath <- combineURI(uri(x), x@genomesFile)
    slots <- slotNames(Genome())
    genomesFields <- c("twoBitPath", "groups", "htmlPath", "metaDb", "trackDb" , "metaTab")
    genomes <- vapply(genomeContainer, function(y) {
        uri <- combineURI(uri(x), y@genome)
        if (!uriExists(uri))
            createResource(uri, dir = TRUE)
        genome <- vapply(slots, function(slotName) {
            slotValue <- slot(y, slotName)
            if (!isEmpty(slotValue) && !is.na(slotValue)) {
                if (slotName %in% genomesFields) {
                    filePath <- combineURI(uri(x), slotValue)
                    if (!uriExists(filePath)) {
                        createResource(filePath)
                    }
                }
                paste0(slotName, " ", slotValue)
            }
            else ""
        }, character(1L))
    }, character(11L))
    genomes <- genomes[genomes != ""]
    genomes <- gsub("\\bgenome\\b", "\ngenome", genomes)
    writeLines(genomes, genomesFilePath)
}

getGenome <- function(x, name) {
    genome <- x@genomeContainer[[name]]
    if (length(genome) == 1L) genome
    else if (length(genome) == 0L) stop("Genome '", name, "' does not exist")
    else if (length(genome) > 1L) stop("Multiple genomes match ", name)
}

setGenome <- function(x, name, value) {
    stopIfNotGenome(value)
    genome <- x@genomeContainer[[name]]
    if (length(genome) == 1L) genome <- value
    else if (length(genome) == 0L) stop("Genome '", name, "' does not exist")
    else if (length(genome) > 1L) stop("Multiple genomes match ", name)
    x@genomeContainer[[name]] <- genome
    x
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
    genomes <- x@genomeContainer
    names(genomes)
})

setMethod("getListElement", "TrackHub", function(x, i, exact = TRUE) {
    TrackHubGenome(x, i)
})

setMethod("names", "TrackHub", function(x) genome(x))

setMethod("length", "TrackHub", function(x) length(names(x)))

setMethod("genomeInfo", "TrackHub", function(x, name) {
    names <- names(x@genomeContainer)
    genome <- x@genomeContainer[[name]]
    if (length(genome) == 0L) stop("Genome '", name, "' does not exist")
    else genome
})

setReplaceMethod("genomeInfo", "TrackHub", function(x, value) {
    stopIfNotGenome(value)
    name <- value@genome
    names <- names(x@genomeContainer)
    genome <- x@genomeContainer[[name]]
    if (length(genome) == 1L) stop("NOTE: Genome '", name, "' already exists")
    else if (length(genome) > 1L) stop("Multiple genomes match ", name)
    else {
        if (!identical(value, Genome())) {
            if (length(x@genomeContainer) == 0L) genomes <- value
            else genomes <- c(unlist(x@genomeContainer), value)
            x@genomeContainer <- GenomeContainer(genomes)
        }
    }
    x
})

setMethod("genomeField", "TrackHub", function(x, name, key) {
    genome <- getGenome(x, name)
    slot(genome, key)
})

setReplaceMethod("genomeField", "TrackHub", function(x, name, key, value) {
    genome <- getGenome(x, name)
    slot(genome, key) <- value
    setGenome(x, name, genome)
})

setMethod("writeTrackHub", "TrackHub", function(x) {
    stopIfNotLocal(hubFile(x))
    setHubContent(x)
    genomesFilePath <- combineURI(uri(x), genomesFile(x))
    if (!uriExists(genomesFilePath) && !is.na(genomesFile(x)))
        createResource(genomesFilePath)
    if (uriExists(genomesFilePath))
        setGenomesContent(x, x@genomeContainer)
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
        }
        else createResource(uri, dir = TRUE)
    }
    th <- new("TrackHub")
    th@uri <- checkURI(uri)
    if (create && !uriExists(hubFile(th))) {
        createResource(hubFile(th))
    } else {
        th <- getHubContent(th)
        genomes <- getGenomesContent(th)
        if (is.list(genomes) && length(genomes) >= 1L)
            th@genomeContainer <- GenomeContainer(unlist(genomes))
    }
    th
}

setAs("character", "TrackHub", function(from) TrackHub(from))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHubGenome class
###

setGeneric("getTracks", function(x) standardGeneric("getTracks"))
setGeneric("trackField", function(x, name, key) standardGeneric("trackField"))
setGeneric("trackField<-", function(x, name, key, value) standardGeneric("trackField<-"))

setClass("TrackHubGenome",
         representation(trackhub = "TrackHub",
                        genome = "character",
                        tracks = "TrackContainer",
                        levels = "integer"),
         contains = "TrackDb")

trackhub <- function(x) x@trackhub

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
    trackDf$value <- gsub("\\b[Oo]n\\b", "TRUE", trackDf$value)
    trackDf$value <- gsub("\\b[Oo]ff\\b", "FALSE", trackDf$value)
    extrafields <- setdiff(trackDf$field, names(fieldToType))
    selectedRows <- !(trackDf$field == extrafields)
    if (length(selectedRows))
        trackDf <- trackDf[selectedRows,]
    args <- Map(as, trackDf$value, fieldToType[trackDf$field])
    names(args) <- trackDf$field
    track <- do.call(Track, args)
    track
}

getTabCountList <- function(filepath) {
    fileContent <- readLines(filepath, warn = FALSE)
    fileContent <- gsub("^\\t*\\s*#(.)*", "", fileContent)
    fileContent <- fileContent[fileContent != ""]
    matches <- gregexpr("^\\s+", fileContent)
    tabCountList <- vapply(matches, attr, integer(1L), "match.length")
    indexes <- grep(-1, tabCountList)
    tabCountList <- replace(tabCountList, indexes, 0)
    as.integer(tabCountList)
}

readFileAsDf <- function(filepath) {
    fileContent <- readLines(filepath, warn = FALSE)
    fileContent <- gsub("^(\\t)*#(.)*", "", fileContent) # to avoid reading commented tracks
    regex <- "^\\s*\\t*([a-zA-Z.]+)\\s?(.*)$"
    field <- sub(regex, "\\1", fileContent)
    value <- sub(regex, "\\2", fileContent)
    contentDf <- data.frame(field, value)
    contentDf <- contentDf[contentDf$field != "",]
}

getTrackDbContent <- function(x, trackDbFilePath) {
    if (isFileEmpty(trackDbFilePath)) {
        x@tracks <- TrackContainer()
        return(x)
    }
    contentDf <- readFileAsDf(trackDbFilePath)
    levels <- getTabCountList(trackDbFilePath)
    tracksIndex <- grep("\\btrack\\b", contentDf$field)
    levels <- levels[tracksIndex]
    totalTracks <- length(tracksIndex)
    tracksIndex[length(tracksIndex) + 1] <- length(contentDf$field) + 1 # to read last track from file
    contentDf$field <- gsub("^(\\t)+", "", contentDf$field)
    # to speed up, reading track by track
    tracks <- lapply(seq_len(totalTracks), function(x) {
        startPosition <- tracksIndex[x]
        endPosition <- tracksIndex[x + 1] - 1
        trackDf <- setNames(data.frame(contentDf$field[startPosition:endPosition],
                                       contentDf$value[startPosition:endPosition]),
                          c("field", "value"))
        track <- createTrack(trackDf)
    })
    x@tracks <- TrackContainer(tracks)
    x@levels <- levels
    x
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

setMethod("trackField", "TrackHubGenome", function(x, name, key) {
    names <- names(x@tracks)
    track <- x@tracks[names == name]
    if (length(track) == 0L) stop("Track '", name, "' does not exist")
    else if (length(track) > 1L) stop("Multiple tracks match ", name)
    slot(track[[1L]], key)
})

setReplaceMethod("trackField", "TrackHubGenome", function(x, name, key, value) {
    names <- names(x@tracks)
    track <- x@tracks[names == name]
    slot(track[[1L]], key) <- value
    x@tracks[names == name] <- track
    x
})

setMethod("organism", "TrackHubGenome", function(object) {
    genome <- getGenome(trackhub(object), genome(object))
    genome@organism
})

setMethod("referenceSequence", "TrackHubGenome", function(x) {
    trackhub <- trackhub(x)
    genome <- getGenome(trackhub, genome(x))
    twoBitPathValue <- genome@twoBitPath
    if (!isFieldEmpty(twoBitPathValue)) {
        twoBitFilePath <- combineURI(uri(trackhub), twoBitPathValue)
        import(twoBitFilePath)
    }
    else stop("genome.txt: 'twoBitPath' does not contain a reference to a file")
})

setReplaceMethod("referenceSequence", "TrackHubGenome", function(x, value) {
    trackhub <- trackhub(x)
    genomesFilePath <- combineURI(uri(trackhub), trackhub@genomesFile)
    stopIfNotLocal(genomesFilePath)
    genome <- getGenome(trackhub, genome(x))
    twoBitPathValue <- genome@twoBitPath
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
    trackhub <- trackhub(x)
    stopIfNotLocal(hubFile(trackhub))
    genome <- getGenome(trackhub, genome(x))
    trackDbValue <- genome@trackDb
    trackDbFilePath <- combineURI(uri(trackhub), trackDbValue)
    if (!isFileEmpty(trackDbFilePath) || length(x@tracks)) {
        tabStrings <- vapply(x@levels, function(y) {
            paste(rep("\t", y), collapse = "")
        },character(1L))
        if (length(tabStrings) == 0)
            tabStrings <- rep("", length(x@tracks))
        slots <- slotNames(Track())
        tracks <- vapply(seq_len(length(x@tracks)), function(i) {
            track <- vapply(slots, function(slotName) {
                slotValue <- slot(x@tracks[[i]], slotName)
                if (!isEmpty(slotValue)) {
                    if (is.na(tabStrings[i])) tabStrings[i] <- ""
                    trackline <- paste0(tabStrings[i], slotName, " ", slotValue)
                    if (slotName == "track") trackline <- paste0("\n", trackline)
                    trackline
                }
                else ""
            }, character(1L))
        }, character(155L))
        tracks <- tracks[tracks != ""]
        tracks <- gsub("\\bTRUE\\b", "on", tracks)
        tracks <- gsub("\\bFALSE\\b", "off", tracks)
        writeLines(tracks, trackDbFilePath)
    }
})

setMethod("show", "TrackHubGenome", function(object) {
    cat(class(object), "track database\ngenome:", genome(object), "\ntrackhub:",
        uri(trackhub(object)), "\n")
    cat(S4Vectors:::labeledLine("names", names(object)))
})

TrackHubGenome <- function(trackhub, genome, create = FALSE) {
    trackhub <- as(trackhub, "TrackHub")
    thg <- new("TrackHubGenome")
    thg@trackhub <- trackhub
    thg@genome <- genome
    genome <- getGenome(trackhub(thg), genome(thg))
    trackDbValue <- genome@trackDb
    if (!isFieldEmpty(trackDbValue)) {
        trackDbFilePath <- combineURI(uri(trackhub(thg)), trackDbValue)
        if (!uriExists(trackDbFilePath) && create) {
            createResource(trackDbFilePath)
        }else if (!isFileEmpty(trackDbFilePath) && uriExists(trackDbFilePath)) {
            thg <- getTrackDbContent(thg, trackDbFilePath)
        }
    }
    thg
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import of tracks from Track Hub
###

setMethod("track", "TrackHubGenome", function(object, name, ...) {
    names <- names(object@tracks)
    track <- object@tracks[names == name]
    if (length(track) == 0L) stop("Track '", name, "' does not exist")
    else if (length(track) > 1L) stop("Multiple tracks match ", name)
    bigDataUrl <- track[[1L]]@bigDataUrl
    parsed <- .parseURI(bigDataUrl)
    if (isEmpty(bigDataUrl)) {
        stop("Track '", name, "' does not contain any data file")
    } else if (uriIsLocal(parsed)) {
        uri <- uri(trackhub(object))
        import(paste0(uri, "/", bigDataUrl), ...)
    } else {
        import(bigDataUrl, ...)
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
    trackhub <- trackhub(object)
    object_uri <- .parseURI(uri(trackhub))
    if (uriIsLocal(object_uri)) {
        genome <- getGenome(trackhub, genome(object))
        trackDbValue <- genome@trackDb
        trackDbValue <- sub(basename(trackDbValue), "", trackDbValue)
        trackDbValue <- sub("/$", "", trackDbValue)
        dest_file <- paste(object_uri$path, trackDbValue, filename, sep = "/")
        dest_file <- sub("^/", "", dest_file)
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
                 signature(object = "TrackHubGenome", value = "BiocFile"),
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
                     genome <- getGenome(trackhub(object), genome(object))
                     trackDbValue <- genome@trackDb
                     trackDbValue <- sub(basename(trackDbValue), "", trackDbValue)
                     trackDbValue <- sub("/$", "", trackDbValue)
                     bigDataUrlValue <- paste(trackDbValue, filename, sep = "/")
                     bigDataUrlValue <- sub("^/", "", bigDataUrlValue)
                     names <- names(object@tracks)
                     trackPosition <- which(names == name)
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

isFileEmpty <- function(path) {
    url <- .parseURI(path)
    if (uriIsLocal(url)) {
        size <- file.size(url$path)
        if (size == 1L || size == 0L)
            TRUE
        FALSE
    } else {
        response <- getURL(path, nobody = T, header = T)
        header <- strsplit(response, "\r\n")[[1]]
        position <- grep("Content-Length:", header)
        contentLength <- sub("Content-Length: ", "", header[position])
        if (contentLength != "NA")
            TRUE
        FALSE
    }
}
