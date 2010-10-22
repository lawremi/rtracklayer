# UCSC genome browser interface

# every UCSC session is identified by a 'hgsid'
setClass("UCSCSession",
         representation(url = "character", hguid = "numeric",
                        views = "environment"),
         contains = "BrowserSession")

# gets an 'hgsid' to initialize the session
setMethod("initialize", "UCSCSession",
          function(.Object, url = "http://genome.ucsc.edu/cgi-bin/",
                   user = NULL, session = NULL, ...)
          {
            .Object@url <- url
            .Object@views <- new.env()
            handle <- getCurlHandle(...)
            gw <- getURL(ucscURL(.Object, "gateway"), cookiefile = tempfile(),
                         header = TRUE, curl = handle)
            cookie <- grep("Set-Cookie: hguid=", gw)
            if (!length(cookie))
              stop("Failed to obtain 'hguid' cookie")
            hguid <- gsub(".*Set-Cookie: hguid=([^;]*);.*", "\\1", gw)
            .Object@hguid <- as.numeric(hguid)
            if (!is.null(user) && !is.null(session)) { ## bring in other session
              ucscGet(.Object, "tracks",
                      list(hgS_doOtherUser = "submit", hgS_otherUserName = user,
                           hgS_otherUserSessionName = session))
            }
            .Object
          })

setMethod("seqlengths", "UCSCSession", function(x) {
  chromInfo <- ucscGet(x, "tracks", list(chromInfoPage = ""))
  path <- "//table/tr[1]/td[contains(text(), 'Sequence name')]/../.."
  table <- getNodeSet(chromInfo, path)[[1]]
  ans <- sapply(getNodeSet(table, "tr/td[@align = 'RIGHT']/text()"), xmlValue)
  ans <- as.integer(gsub("[^0-9]", "", ans))
  names(ans) <- sapply(getNodeSet(table, "tr/td/a/text()"), xmlValue)
  ans
})

setMethod("seqnames", "UCSCSession", function(x) names(seqlengths(x)))

normArgTrackData <- function(value, session) {
  genomes <- lapply(value, genome)
  genomes[sapply(genomes, length) == 0L] <- ""
  tapply(value, unlist(genomes),
         function(tracks)
         {
           genome <- genome(tracks[[1]])
           if (length(genome))
             genome(session) <- genome
           spaces <- unlist(lapply(tracks, names))
           badSpaces <- setdiff(spaces, seqnames(session))
           if (length(badSpaces))
             stop("Invalid chromosomes for ", genome(session), ": ",
                  paste(badSpaces, collapse = ", "))
         })
  value
}

setReplaceMethod("track", c("UCSCSession", "RangedDataList"),
          function(object, name = names(value),
                   format = c("auto", "bed", "wig", "gff1", "bed15"), ...,
                   value)
          {
            format <- match.arg(format)
            if (length(value)) {
              ## upload values in blocks, one for each genome
              value <- normArgTrackData(value, session)
              names(value) <- name
              genomes <- lapply(value, genome)
              genomes[sapply(genomes, length) == 0L] <- ""
              tapply(value, unlist(genomes),
                     function(tracks)
                     {
                       form <- ucscForm(tracks, format, ...)
                       response <- ucscPost(object, "custom", form)
### FIXME: need to check for error
                     })
            }
            object
          })

setMethod("browserViews", "UCSCSession",
          function(object) object@views$instances)

## get the list of track names
setMethod("trackNames", "UCSCSession",
          function(object) ucscTracks(object)@ids)

## get the current range
setMethod("range", "UCSCSession",
          function(x, ..., na.rm) range(ucscCart(x)))

setReplaceMethod("range", "UCSCSession",
                 function(x, value) {
                   ucscGet(x, "cart", ucscForm(normGenomeRange(value, x)))
                   x
                 })

setReplaceMethod("genome", "UCSCSession",
                 function(x, value) {
                   ucscGet(x, "gateway", list(db = value))
                   x
                 })

GRangesForUCSCGenome <- function(genome, chrom = NULL, ranges = NULL, ...)
{
  if (missing(genome) || !IRanges:::isSingleString(genome))
    stop("'genome' must be a single string identifying a genome")
  session <- browserSession("UCSC")
  genome(session) <- genome
  GRangesForGenome(genome, seqlengths(session), chrom = chrom, ranges = ranges,
                   ...)
}


## context for querying UCSC tables
setClass("UCSCTableQuery",
         representation(session = "UCSCSession",
                        track = "characterORNULL",
                        table = "characterORNULL",
                        range = "GRanges",
                        outputType = "characterORNULL",
                        NAMES = "characterORNULL",
                        intersectTrack = "characterORNULL"))

setMethod("show", "UCSCTableQuery",
          function(object) {
            cat("Get ")
            if (!is.null(tableName(object)))
              cat("table '", tableName(object), "' from ", sep = "")
            cat("track '", names(trackName(object)), "' within ", sep = "")
            range <- range(object)
            if (length(range) > 1)
              start <- end <- chrom <- "*"
            else {
              chrom <- as.character(seqnames(range))
              start <- start(range)
              end <- end(range)
            }
            cat(genome(range), ":", chrom, ":", start, "-", end, sep="")
            for (itrack in names(intersectTrack(object)))
              cat(" &", itrack)
            cat("\n")
          })

setMethod("browserSession", "UCSCTableQuery", function(object) {
  object@session
})

setGeneric("browserSession<-",
           function(object, ..., value) standardGeneric("browserSession<-"))
setReplaceMethod("browserSession", c("UCSCTableQuery", "UCSCSession"),
                 function(object, value) {
                   object@session <- value
                   object
                 })

setMethod("range", "UCSCTableQuery", function(x, ..., na.rm) x@range)
setReplaceMethod("range", "UCSCTableQuery",
                 function(x, value) {
                   x@range <- normGenomeRange(value, browserSession(x))
                   x
                 })

setGeneric("trackName", function(x, ...) standardGeneric("trackName"))
setMethod("trackName", "UCSCTableQuery", function(x) x@track)

setGeneric("trackName<-",
           function(x, ..., value) standardGeneric("trackName<-"))
setReplaceMethod("trackName", "UCSCTableQuery", function(x, value)
                 {
                   x@track <- normArgTrack(value, x)
                   x
                 })

setGeneric("tableName", function(x, ...) standardGeneric("tableName"))
setMethod("tableName", "UCSCTableQuery", function(x) x@table)

normArgTable <- function(name, query) {
  if (!is.null(name)) {
    if (!isSingleString(name))
      stop("table name must be a single string or NULL")
    if (!name %in% tableNames(query))
      stop("unknown table name '", name, "'")
  }
  name
}

setGeneric("tableName<-", function(x, ..., value)
           standardGeneric("tableName<-"))
setReplaceMethod("tableName", "UCSCTableQuery", function(x, value)
                 {
                   x@table <- normArgTable(value, x)
                   x
                 })

setMethod("names", "UCSCTableQuery", function(x) x@NAMES)
setReplaceMethod("names", "UCSCTableQuery", function(x, value) {
  x@NAMES <- value
  x
})

setGeneric("intersectTrack", function(x, ...)
           standardGeneric("intersectTrack"))
setMethod("intersectTrack", "UCSCTableQuery", function(x) x@intersectTrack)
setGeneric("intersectTrack<-", function(x, ..., value)
           standardGeneric("intersectTrack<-"))
setReplaceMethod("intersectTrack", "UCSCTableQuery", function(x, value) {
  x@intersectTrack <- normArgTrack(value, x, TRUE)
  x
})

## not exported
setGeneric("outputType", function(x, ...) standardGeneric("outputType"))
setMethod("outputType", "UCSCTableQuery", function(x) x@outputType)
setGeneric("outputType<-",
           function(x, ..., value) standardGeneric("outputType<-"))
setReplaceMethod("outputType", "UCSCTableQuery",
                 function(x, value) {
                   x@outputType <- value
                   x
                 })

normArgTrack <- function(name, trackids) {
  if (is.null(name))
    return(name)
  if (!isSingleString(name))
    stop("'track' must be a single string")
  if (is(trackids, "UCSCTableQuery"))
    trackids <- trackNames(trackids)
  if (!(name %in% trackids)) {
    mapped_name <- trackids[name]
    if (is.na(mapped_name))
      stop("Unknown track: ", name)
    name <- mapped_name
  } else names(name) <- name
  name
}

setGeneric("ucscTableQuery", function(x, ...) standardGeneric("ucscTableQuery"))
setMethod("ucscTableQuery", "UCSCSession",
          function(x, track = NULL, range = genome(x), table = NULL,
                   names = NULL, intersectTrack = NULL)
          {
            if (!is(names, "characterORNULL"))
              stop("'names' must be 'NULL' or a character vector")
            ## only inherit the genome from the session
            range <- normGenomeRange(range, x)
            query <- new("UCSCTableQuery", session = x, range = range,
                         NAMES = names)
            ## the following line must always happen to initialize the session
            ## otherwise stuff can go haywire
            trackids <- trackNames(query)
            if (!is.null(track) || !is.null(intersectTrack)) {
              query@track <- normArgTrack(track, trackids)
              query@intersectTrack <- normArgTrack(intersectTrack, trackids)
            }
            tableName(query) <- table
            query
          })

ucscTableGet <- function(query, .parse = TRUE, tracks = FALSE, ...)
  ucscGet(browserSession(query), "tables",
          c(ucscForm(query, tracks = tracks), ...), .parse = .parse)

## gets the track names available from the table browser

setMethod("trackNames", "UCSCTableQuery",
          function(object) {
            doc <- ucscTableGet(object, tracks = TRUE)
            track_path <- "//select[@name = 'hgta_track']/option/@value"
            tracks <- unlist(getNodeSet(doc, track_path))
            label_path <- "//select[@name = 'hgta_track']/option/text()"
            labels <- sub("\n.*$", "",
                          sapply(getNodeSet(doc, label_path), xmlValue))
            names(tracks) <- labels
            tracks
          })

## returns a character vector of table names for a given track name + range
setGeneric("tableNames", function(object, ...)
           standardGeneric("tableNames"))

setMethod("tableNames", "UCSCTableQuery",
          function(object, trackOnly = FALSE)
          {
            doc <- ucscTableGet(object)
            table_path <- "//select[@name = 'hgta_table']/option/@value"
            tables <- unlist(getNodeSet(doc, table_path))
            outputType <- outputType(object)
            if (trackOnly) {
              trackOutputs <- c("wigData", "wigBed", "bed")
              if (!is.null(outputType))
                outputType <- intersect(trackOutputs, outputType)
              else outputType <- trackOutputs
            }
            if (!is.null(outputType)) {
              checkOutput <- function(table) {
                object@table <- table # avoid accessor to skip check
                outputs <- ucscTableOutputs(object)
                any(outputType %in% outputs)
              }
              tables <- tables[sapply(tables, checkOutput)]
            }
            unname(tables)
          })

setGeneric("ucscTableOutputs",
           function(object, ...)
           standardGeneric("ucscTableOutputs"))

## returns a character vector of available output types for the table
## not exported
setMethod("ucscTableOutputs", "UCSCTableQuery",
          function(object) {
            doc <- ucscTableGet(object)
            output_path <- "//select[@name = 'hgta_outputType']/option/@value"
            unlist(getNodeSet(doc, output_path))
          })

setGeneric("ucscTableSchema",
           function(object, ...) standardGeneric("ucscTableSchema"))

setMethod("ucscTableSchema", "UCSCTableQuery",
          function(object)
          {
            doc <- ucscTableGet(object, hgta_doSchema = "describe table schema")
### TODO: get the schema as a data.frame for the given track and table 
          })

## export data from UCSC (internal utility)
ucscExport <- function(object)
{
  get_hgsid <- function(node)
    getNodeSet(node, "//input[@name = 'hgsid']/@value")[[1]]
  hgsid <- NULL
  if (!is.null(names(object))) { # filter by names
    text <- paste(names(object), collapse = "\n")
    output <- ucscTableGet(object, hgta_doPastedIdentiers = "submit",
                           hgta_pastedIdentifiers = text)
    error <- getNodeSet(output,
                        "//script[contains(text(), '{showWarnBox')]/text()")
    if (length(error))
      warning(sub(".*'<li>(.*?)'.*", "\\1", xmlValue(error[[1]])))
    hgsid <- get_hgsid(output)
  }
  if (!is.null(intersectTrack(object))) {
    itrack <- intersectTrack(object)
    iquery <- object
    iquery@track <- itrack
    itable <- tableNames(iquery, TRUE)
    if (!length(itable))
      stop("No table for intersection track: ", itrack)
    if (length(itable) > 1) # for now anyway
      itable <- itable[1]
    output <- ucscTableGet(object, hgta_nextIntersectGroup = "allTracks",
                           hgta_nextIntersectTrack = itrack,
                           hgta_nextIntersectTable = itable,
                           hgta_nextIntersectOp = "any",
                           hgta_doIntersectSubmit = "submit",
                           boolshad.hgta_nextInvertTable = "0",
                           boolshad.hgta_nextInvertTable2 = "0",
                           hgsid = hgsid)
    hgsid <- get_hgsid(output)
  }
  followup <- NULL
  if (outputType(object) == "bed") { ## some formats have extra pages
    followup <- list(hgta_doGetBed = "get BED",
                     hgta_printCustomTrackHeaders = "on",
                     boolshad.hgta_printCustomTrackHeaders = "1")
  }
  output <- ucscTableGet(object, !is.null(followup),
                         hgta_doTopSubmit = "get output",
                         hgsid = hgsid)
  if (!is.null(followup)) {
    hgsid <- get_hgsid(output)
    form <- c(followup, list(hgsid = hgsid))
    output <- ucscGet(browserSession(object), "tables", form, .parse = FALSE)
  }
  output
}

setMethod("track", "UCSCSession",
          function(object, name, range = base::range(object), table = NULL,
                   asRangedData = TRUE)
          {
            track(ucscTableQuery(object, name, range, table), asRangedData)
          })

## download a trackSet by name
setMethod("track", "UCSCTableQuery",
          function(object, asRangedData = TRUE)
          {
            tables <- tableNames(object)
            table <- tableName(object)
            if (!is.null(table) && !(table %in% tables))
              stop("Unknown table: '", table, "'. Valid table names: ", tables)
            formats <- c("wigData", "wigBed", "bed")
            ## attempt to automatically determine the table
            if (!is.null(table))
              tables <- table
            for (table in tables) {
              outputs <- ucscTableOutputs(object)
              if (any(formats %in% outputs))
                break
            }
            if (!any(formats %in% outputs))
              stop("No supported output types")
            if ("wigData" %in% outputs) { # track stored as wig
              format <- "wig"
              output <- "wigData"
            } else {
              format <- output <- "bed"
              if ("wigBed" %in% outputs)
                output <- "wigBed"
            }
            outputType(object) <- output
            tableName(object) <- table
            output <- ucscExport(object)
            import(text = output, format = format, asRangedData = asRangedData)
          })

## grab sequences for features in 'track' at 'range'
## setMethod("getSeq", "UCSCSession",
##           function(object, range, table = "gold")
##           {
##             followup <- list(hgta_doGenomicDna = "get sequence",
##                              hgSeq.casing = "upper",
##                              hgSeq.repMasking = "lower")
##             output <- ucscExport(object, range, "gold", table, "sequence",
##                                  followup)
##             con <- file()
##             writeLines(output, con)
##             set <- read.DNAStringSet(con, "fasta")
##             close(con)
##             set
##           })

## get a data.frame from a UCSC table
## think about taking specific columns
setGeneric("getTable",
           function(object, ...) standardGeneric("getTable"))
setMethod("getTable", "UCSCTableQuery",
          function(object)
          {
            if (!("primaryTable" %in% ucscTableOutputs(object)))
              stop("tabular output format not available")
            outputType(object) <- "primaryTable"
            if (is.null(tableName(object))) # must specify a table name
              tableName(object) <- tableNames(object)[1]
            output <- ucscExport(object)
            ## since '#' is not treated as a comment, we discard the
            ## error message, leaving only the header
            if (grepl("\\n# No results", output))
              output <- gsub("\\n.*", "", output)
            f <- file()
            writeLines(output, f)
            header <- readChar(f, 1) ## strip off the '#' header prefix
            tab <- read.table(f, sep = "\t", header=TRUE, comment.char = "",
                              quote = "")
            close(f)
            tab
          })

## UCSC genome view
setClass("UCSCView", representation(hgsid = "numeric"),
         contains = "BrowserView")

## create a view for the given session, position and track visibility settings
## if 'tracks' is a character vector (but not a UCSCTrackModes instance) it is
## assumed to name the tracks that should be in the view. otherwise, an
## attempt is made to coerce it to a UCSCTrackModes instance.
setMethod("browserView", "UCSCSession",
          function(object, range, track, imagewidth = 800, ...)
          {
            form <- list()
            if (!missing(range)) {
              if (length(range) > 1) {
                ranges <- range
                views <- vector("list", length(ranges))
                for (i in seq(length(ranges))) {
                  range <- ranges[i]
                  views[[i]] <- callGeneric()
                }
                return(BrowserViewList(views))
              }
              range <- normGenomeRange(range, object)
              form <- c(form, ucscForm(range))
            }
            view <- new("UCSCView", session = object)
            ## new hgsid for each browser launch
            doc <- ucscGet(object, "gateway")
            node <- getNodeSet(doc, "//input[@name = 'hgsid']/@value")[[1]]
            hgsid <- node ##xmlValue(node)
            view@hgsid <- as.numeric(hgsid)
            ## figure out track modes
            origModes <- modes <- ucscTrackModes(view)
            if (!missing(track)) {
              if (class(track) == "character")
                trackNames(modes) <- track
              else {
                userModes <- as(track, "UCSCTrackModes")
                modes[names(userModes)] <- userModes
              }
            }
            argModes <- ucscTrackModes(...)
            modes[names(argModes)] <- argModes
            modes <- modes[modes != origModes]
            form <- c(form, ucscForm(modes), ucscForm(view))
            if (!missing(imagewidth))
              form <- c(form, pix = imagewidth)
            ## launch a web browser
            ucscShow(object, "tracks", form)
            ## save this view
            object@views$instances <- c(object@views$instances, view)
            view
          })

# every view has a "mode" (hide, dense, pack, squish, full) for each track
### FIXME: probably should be merged with ucscTracks
### Or just leave it; ucscTracks might become more complex, while we
### need a simple way to manipulate track modes.
setClass("UCSCTrackModes", representation(labels = "character"),
         contains = "character")

# get/set track modes to/from e.g. a view
setGeneric("ucscTrackModes",
           function(object, ...) standardGeneric("ucscTrackModes"))

# convenience constructor for track mode object
setMethod("ucscTrackModes", "character",
          function(object, labels, hide = character(),
                   dense = character(), pack = character(),
                   squish = character(), full = character())
          {
            object[hide] <- "hide"
            object[dense] <- "dense"
            object[pack] <- "pack"
            object[squish] <- "squish"
            object[full] <- "full"
            if (missing(labels))
              labels <- names(object)
            new("UCSCTrackModes", object, labels = as.character(labels))
          })
setMethod("ucscTrackModes", "missing",
          function(object, ...) ucscTrackModes(character(), ...))

setMethod("ucscTrackModes", "UCSCView",
          function(object)
          {
            ucscTrackModes(ucscTracks(object))
          })

setMethod("ucscTrackModes", "UCSCSession",
          function(object)
          {
            ucscTrackModes(ucscTracks(object))
          })

setGeneric("ucscTrackModes<-",
           function(object, value) standardGeneric("ucscTrackModes<-"))
setReplaceMethod("ucscTrackModes", c("UCSCView", "UCSCTrackModes"),
                 function(object, value)
                 { # FIXME: needs to be more extensible
                   browserView(object@session, range(object), value)
                 })
setReplaceMethod("ucscTrackModes", c("UCSCView", "character"),
                 function(object, value)
                 {
                   ucscTrackModes(object) <- ucscTrackModes(value)
                   object
                 })

## subsetting UCSCTrackModes

## if not in ids, try labels
resolveTrackIndex <- function(object, i) {
  if (is.character(i)) {
    missing <- !(i %in% names(object))
    matching <- match(i[missing], object@labels)
    if (any(is.na(matching))) {
      unmatched <- i[missing][is.na(matching)]
      stop("Unknown track(s): ", paste(unmatched, collapse = ", "))
    }
    i[missing] <- names(object)[matching]
  }
  i
}

setMethod("[", "UCSCTrackModes", function(x, i, j, ..., drop=FALSE) {
  vec <- x@.Data
  names(vec) <- names(x)
  names(x@labels) <- names(x)
  ind <- resolveTrackIndex(x, i)
  initialize(x, vec[ind], labels = x@labels[ind])
})

setReplaceMethod("[", "UCSCTrackModes", function(x, i, j, ..., value) {
  vec <- x@.Data
  names(vec) <- names(x)
  vec[resolveTrackIndex(x, i)] <- value
  x@.Data <- as.vector(vec)
  x
})

# handle simple track show/hide

setMethod("trackNames", "UCSCTrackModes",
          function(object)
          {
            visible <- object != "hide"
            tracks <- names(object)[visible]
            names(tracks) <- object@labels[visible]
            tracks
          })
setReplaceMethod("trackNames", "UCSCTrackModes",
                 function(object, value)
                 {
                   value <- resolveTrackIndex(object, value)
                   spec <- names(object) %in% value
                   object[!spec] <- "hide"
                   object[spec & object == "hide"] <- "full"
                   object
                 })

setMethod("trackNames", "UCSCView",
          function(object)
          {
            tracks <- ucscTracks(object)
            modes <- ucscTrackModes(tracks)
            tracks@ids[tracks@ids %in% trackNames(modes)]
          })
setReplaceMethod("trackNames", "UCSCView",
                 function(object, value)
                 {
                   trackNames(ucscTrackModes(object)) <- value
                   object
                 })


setMethod("visible", "UCSCView", function(object) {
  modes <- ucscTrackModes(object)
  vis <- modes != "hide"
  names(vis) <- modes@labels
  vis
})
setReplaceMethod("visible", "UCSCView", function(object, value) {
  modes <- ucscTrackModes(object)
  modes[value & modes == "hide"] <- "full"
  modes[!value] <- "hide"
  ucscTrackModes(object) <- modes
  object
})


setMethod("range", "UCSCView",
          function(x, ..., na.rm) range(ucscCart(x)))
setReplaceMethod("range", "UCSCView",
                 function(x, value)
                 {
                   browserView(x@session, value, ucscTrackModes(x))
                 })

# only one view per session; a view is always active
setMethod("activeView", "UCSCView", function(object) TRUE)

# ucscTrackSet

# visual properties are specified by a "track line" for UCSC
setClass("TrackLine",
         representation(name = "character", description = "character",
                        visibility = "character", color = "integer",
                        priority = "numeric"),
         prototype(name = "R Track"))

setMethod("show", "TrackLine",
          function(object)
          {
            cat(as(object, "character"), "\n")
          })

setClass("BasicTrackLine",
         representation(itemRgb = "logical", useScore = "logical",
                        group = "character", db = "character",
                        offset = "numeric", url = "character",
                        htmlUrl = "character"),
         contains = "TrackLine")

ucscPair <- function(key, value) paste(key, value, sep = "=")

# to a text line
setAs("TrackLine", "character",
      function(from)
      {
        checkString <- function(str, len) {
          ## These are more annoying than useful
          ## if (nchar(gsub("[a-zA-Z0-9_ ]", "", str)))
          ##   warning("The string '", str,
          ##           "' contains non-standard characters.")
          ## if (nchar(str) > len) {
          ##   str <- substring(str, 1, len)
          ##   warning("The string '", str, "' must be less than ", len,
          ##           " characters; it has been truncated.")
          ## }
          if (regexpr(" ", str)[1] != -1)
            str <- paste("\"", str, "\"", sep="")
          str
        }
        str <- "track"
        name <- from@name
        if (length(name))
          str <- paste(str, " name=", checkString(name, 15), sep="")
        desc <- from@description
        if (length(desc))
          str <- paste(str, " description=", checkString(desc, 60), sep="")
        vis <- from@visibility
        if (length(vis))
          str <- paste(str, " visibility=", vis, sep="")
        color <- from@color
        if (length(color))
          str <- paste(str, " color=", paste(color, collapse=","), sep="")
        priority <- from@priority
        if (length(priority))
          str <- paste(str, " priority=", priority, sep="")
        str
      })

setAs("BasicTrackLine", "character",
      function(from)
      {
        str <- as(as(from, "TrackLine"), "character")
        itemRgb <- from@itemRgb
        if (length(itemRgb) && itemRgb)
          str <- paste(str, "itemRgb=On")
        useScore <- from@useScore
        if (length(useScore) && useScore)
          str <- paste(str, "useScore=1")
        group <- from@group
        if (length(group))
          str <- paste(str, " group=", group, sep="")
        db <- from@db
        if (length(db))
          str <- paste(str, " db=", db, sep="")
        offset <- from@offset
        if (length(offset))
          str <- paste(str, " offset=", offset, sep="")
        url <- from@url
        if (length(url))
          str <- paste(str, " url=", url, sep="")
        htmlUrl <- from@htmlUrl
        if (length(htmlUrl))
          str <- paste(str, " htmlUrl=", htmlUrl, sep="")
        str
      })

ucscParsePairs <- function(str)
{
  str <- sub("^[^ ]* ", "", str)
  split <- strsplit(str, "=")[[1]]
  split <- sub("^ +", "", sub(" +$", "", split))
  vals <- character(0)
  if (length(split)) {
    mixed <- tail(head(split, -1), -1)
    tags <- head(split, 1)
    vals <- tail(split, 1)
    if (length(mixed)) {
      tags <- c(tags, sub(".* ([^ ]*)$", "\\1", mixed))
      vals <- c(sub("(.*) [^ ]*$", "\\1", mixed), vals)
    }
    names(vals) <- tags
    vals <- sub("\"?([^\"]*)\"?", "\\1", vals)
  }
  vals
}

# from a text line
setAs("character", "TrackLine",
      function(from)
      {
        line <- new("TrackLine")
        vals <- ucscParsePairs(from)
        if (!is.na(vals["name"]))
          line@name <- vals["name"]
        if (!is.na(vals["description"]))
          line@description <- vals["description"]
        if (!is.na(vals["visibility"]))
          line@visibility <- vals["visibility"]
        if (!is.na(vals["color"]))
          line@color <- as.integer(strsplit(vals["color"], ",")[[1]])
        if (!is.na(vals["priority"]))
          line@priority <- as.numeric(vals["priority"])
        line
      })

setAs("character", "BasicTrackLine",
      function(from)
      {
        line <- new("BasicTrackLine", as(from, "TrackLine"))
        vals <- ucscParsePairs(from)
        if (!is.na(vals["itemRgb"]))
          line@itemRgb <- vals["itemRgb"] == "On"
        if (!is.na(vals["useScore"]))
          line@useScore <- vals["useScore"] == "1"
        if (!is.na(vals["group"]))
          line@group <- vals["group"]
        if (!is.na(vals["db"]))
          line@db <- vals["db"]
        if (!is.na(vals["offset"]))
          line@offset <- vals["offset"]
        if (!is.na(vals["url"]))
          line@url <- vals["url"]
        if (!is.na(vals["htmlUrl"]))
          line@htmlUrl <- vals["htmlUrl"]
        line
      })


setClass("GraphTrackLine",
         representation(altColor = "integer", autoScale = "logical",
                        gridDefault = "logical", maxHeightPixels = "numeric",
                        graphType = "character", viewLimits = "numeric",
                        yLineMark = "numeric", yLineOnOff = "logical",
                        windowingFunction = "character",
                        smoothingWindow = "numeric", type = "character"),
         contains = "TrackLine")

setAs("GraphTrackLine", "character",
      function(from)
      {
        str <- as(as(from, "TrackLine"), "character")
        type <- if (from@type == "wig") "wiggle_0" else "bedGraph"
        str <- paste(str, " type=", type, sep = "")
        color <- from@altColor
        if (length(color))
          str <- paste(str, " altColor=", paste(color, collapse=","), sep="")
        autoScale <- from@autoScale
        if (length(autoScale) && !autoScale)
          str <- paste(str, "autoScale=off")
        gridDefault <- from@gridDefault
        if (length(gridDefault) && gridDefault)
          str <- paste(str, "gridDefault=On")
        maxHeightPixels <- from@maxHeightPixels
        if (length(maxHeightPixels) && maxHeightPixels)
          str <- paste(str, " maxHeightPixels=",
                       paste(maxHeightPixels, collapse=":"), sep = "")
        graphType <- from@graphType
        if (length(graphType))
          str <- paste(str, " graphType=", graphType, sep = "")
        viewLimits <- from@viewLimits
        if (length(viewLimits))
          str <- paste(str, " viewLimits=", paste(viewLimits, collapse = ":"),
                       sep = "")
        yLineMark <- from@yLineMark
        if (length(yLineMark))
          str <- paste(str, " yLineMark=", yLineMark, sep = "")
        yLineOnOff <- from@yLineOnOff
        if (length(yLineOnOff) && yLineOnOff)
          str <- paste(str, "yLineOnOff=On")
        windowingFunction <- from@windowingFunction
        if (length(windowingFunction))
          str <- paste(str, " windowingFunction=", windowingFunction, sep = "")
        smoothingWindow <- from@smoothingWindow
        if (length(smoothingWindow))
          str <- paste(str, " smoothingWindow=", smoothingWindow, sep = "")
        str
      })

setAs("character", "GraphTrackLine",
      function(from)
      {
        line <- new("GraphTrackLine", as(from, "TrackLine"))
        vals <- ucscParsePairs(from)
        type <- vals["type"]
        if (!(type %in% c("wiggle_0", "bedGraph")))
          stop("Unknown graph track type: ", type)
        line@type <- if (type == "wiggle_0") "wig" else "bedGraph"
        if (!is.na(vals["altColor"]))
          line@altColor <- as.integer(strsplit(vals["altColor"], ",")[[1]])
        if (!is.na(vals["autoScale"]))
          line@autoScale <- vals["autoScale"] == "On"
        if (!is.na(vals["gridDefault"]))
          line@gridDefault <- vals["gridDefault"] == "On"
        if (!is.na(vals["maxHeightPixels"]))
          line@maxHeightPixels <-
            as.numeric(strsplit(vals["maxHeightPixels"], ":")[[1]])
        if (!is.na(vals["graphType"]))
          line@graphType <- vals["graphType"]
        if (!is.na(vals["viewLimits"]))
          line@viewLimits <-
            as.numeric(strsplit(vals["viewLimits"], ":")[[1]])
        if (!is.na(vals["yLineMark"]))
          line@yLineMark <- as.numeric(vals["yLineMark"])
        if (!is.na(vals["yLineOnOff"]))
          line@yLineOnOff <- vals["yLineOnOff"] == "On"
        if (!is.na(vals["windowingFunction"]))
          line@windowingFunction <- vals["windowingFunction"]
        if (!is.na(vals["smoothingWindow"]))
          line@smoothingWindow <- as.numeric(vals["smoothingWindow"])
        line
      })

setAs("BasicTrackLine", "GraphTrackLine",
      function(from) new("GraphTrackLine", from))

setAs("GraphTrackLine", "BasicTrackLine",
      function(from) new("BasicTrackLine", from))

setClass("UCSCData",
         representation(trackLine = "TrackLine"),
         prototype(trackLine = new("BasicTrackLine")),
         "RangedData")

setMethod("show", "UCSCData",
          function(object)
          {
            if (!is.null(object@trackLine@name))
              cat("UCSC track '", object@trackLine@name, "'\n", sep = "")
            callNextMethod()
          })

ucscNormSeqNames <- function(nms) {
  nms <- gsub("^([0-9A-Z]+)$", "chr\\1", nms)
  nms
}

chooseGraphType <- function(from) {
  r <- ranges(from)[[1]] # heuristic only needs first chromosome
  type <- "bedGraph"
  ## decide whether compression is a good idea
  steps <- diff(sort(start(r)))
  if (length(unique(steps)) == 1 || # fixed-step makes sense
      ((3 * length(unique(width(r)))) < length(r) && # makes sense wrt size
       mean(steps) < 100)) # dense enough for UCSC efficiency
    type <- "wig"
  type
}

setAs("RangedData", "UCSCData", function(from) {
  names(from) <- ucscNormSeqNames(names(from))
  if (is.numeric(score(from))) { # have numbers, let's plot them
    type <- chooseGraphType(from)
    line <- new("GraphTrackLine", type = type)
  } else line <- new("BasicTrackLine")
  new("UCSCData", from, trackLine = line)
})

# the 'ucsc' format is a meta format with a track line followed by
# tracks formatted as 'wig', 'bed', 'bed15', 'bedGraph', 'gff', 'gtf', or 'psl'.
# 'gtf' and 'psl' are not yet supported
setGeneric("export.ucsc",
           function(object, con,
                    subformat = c("auto", "gff1", "wig", "bed", "bed15",
                      "bedGraph"),
                    append = FALSE, ...)
           standardGeneric("export.ucsc"))

setMethod("export.ucsc", "RangedDataList",
          function(object, con,
                   subformat = c("auto", "gff1", "wig", "bed", "bed15",
                     "bedGraph"),
                   append, trackNames, ...)
          {
            subformat <- match.arg(subformat)
            if (missing(trackNames)) {
              trackNames <- names(object)
              if (is.null(trackNames))
                trackNames <- paste("R Track", seq_len(length(object)))
              ucsc <- unlist(lapply(object, is, "UCSCData"))
              lines <- unlist(lapply(object[ucsc], slot, "trackLine"))
              trackNames[ucsc] <- as.character(sapply(lines, slot, "name"))
            }
            for (i in seq_len(length(object))) {
              export.ucsc(object[[i]], con, subformat, name = trackNames[i],
                          append = append, ...)
              append <- TRUE
            }
          })

trackLineClass <- function(subformat)
{
  if (subformat == "wig" || subformat == "bedGraph")
    "GraphTrackLine"
  else if (subformat == "bed15")
    "Bed15TrackLine"
  else "BasicTrackLine"
}

setMethod("export.ucsc", "ANY",
          function(object, con, subformat, append, ...)
          {
            cl <- class(object)
            track <- try(as(object, "RangedData"), silent = TRUE)
            if (class(track) == "try-error") {
              track <- try(as(object, "RangedDataList"), silent = TRUE)
              if (class(track) == "try-error")
                stop("cannot export object of class '", cl, "'")
            }
            export.ucsc(track, con = con, subformat = subformat,
                        append = append, ...)
          })

setMethod("export.ucsc", "RangedData",
          function(object, con, subformat, append, ...)
          {
            object <- as(object, "UCSCData")
            subformat <- match.arg(subformat)
            export.ucsc(object, con, subformat, append, ...)
           })

setMethod("export.ucsc", c("UCSCData", "characterORconnection"),
          function(object, con, subformat, append, ...)
          {
            subformat <- match.arg(subformat)
            auto <- FALSE
            if (subformat == "auto") {
              auto <- TRUE
              subformat <- "bed"
              if (is(object@trackLine, "Bed15TrackLine"))
                subformat <- "bed15"
              else if (is(object@trackLine, "GraphTrackLine"))
                subformat <- object@trackLine@type
            }
            graphFormat <- subformat %in% c("wig", "bedGraph")
            if (graphFormat) {
              strand <- as.character(strand(object))
              strand[is.na(strand)] <- "NA"
              isStrandDisjoint <- function(track) {
                all(unlist(lapply(ranges(track), function(r) {
                  isDisjoint(r) && all(width(r) > 0)
                })))
              }
              if (!all(unlist(lapply(split(object, strand), isStrandDisjoint))))
              {
                if (auto)
                  subformat <- "bed"
                else stop("Track not compatible with WIG/bedGraph: ",
                          "Overlapping features must be on separate strands",
                          " and every feature width must be positive")
              }
            }
            lineClass <- trackLineClass(subformat)
            if (!is(object@trackLine, lineClass))
              object@trackLine <- as(object@trackLine, lineClass)
            if (is(object@trackLine, "GraphTrackLine"))
              object@trackLine@type <- subformat
            args <- list(...)
            lineArgs <- names(args) %in% slotNames(lineClass)
            for (argName in names(args)[lineArgs])
              slot(object@trackLine, argName) <- args[[argName]]
            trackLine <- NULL
            if (graphFormat) {
              subformatOrig <- subformat
              subformat <- paste(subformat, "Lines", sep = "")
              strand <- as.character(strand(object))
              strand[is.na(strand)] <- "NA"
              if (!all(strand[1] == strand)) {
                nameMap <- c("+" = "p", "-" = "m", "NA" = "NA")
                strand <- factor(strand)
                name <- paste(object@trackLine@name, nameMap[levels(strand)])
                tracks <- split(object, strand)
                export.ucsc(tracks, con, subformatOrig, append,
                            trackNames = name, ...)
                return()
              }
            } else if (subformat == "bed15") {
              subformat <- "bed15Lines"
              if (is.null(object@trackLine@expNames))
                object@trackLine@expNames <- colnames(object)
              trackLine <- object@trackLine
            }
            cat(as(object@trackLine, "character"), "\n", file=con, sep = "",
                append = append)
            do.call(export, c(list(as(object, "RangedData"), con, subformat),
                                args[!lineArgs], trackLine = trackLine,
                              append = TRUE))
          })

# for GFF, the track line should go in a comment
setMethod("export.gff", c("UCSCData", "characterORconnection"),
          function(object, con, version, source, append)
          {
            gffComment(con, as(object@trackLine, "character"))
            callNextMethod()
          })

setGeneric("import.ucsc",
           function(con, subformat =
                   c("auto", "gff1", "wig", "bed", "bed15", "bedGraph"),
                    drop = FALSE, asRangedData = TRUE, ...)
           standardGeneric("import.ucsc"))
setMethod("import.ucsc", "characterORconnection",
          function(con, subformat, drop = FALSE, asRangedData = TRUE, ...)
          {
            if (!IRanges:::isTRUEorFALSE(asRangedData))
              stop("'asRangedData' must be TRUE or FALSE")
            subformat <- match.arg(subformat)
            lines <- readLines(con, warn = FALSE)
            tracks <- grep("^track", lines)
            trackLines <- lines[tracks]
            starts <- tracks + 1L
            ends <- c(tail(tracks, -1) - 1L, length(lines))
            makeTrackSet <- function(i)
            {
              line <- as(trackLines[i], trackLineClass(subformat))
              if (subformat == "wig" || subformat == "bedGraph")
                subformat <- paste(subformat, "Lines", sep = "")
              if (starts[i] <= ends[i])
                text <- window(lines, starts[i], ends[i])
              else
                text <- character()
              if (subformat == "bed15") # need to pass track line
                ucsc <- import(format = "bed15Lines", text = text,
                               trackLine = line,
                               asRangedData = asRangedData, ...)
              else
                ucsc <- import(format = subformat, text = text,
                               asRangedData = asRangedData, ...)
              if (asRangedData) {
                ucsc <- as(ucsc, "UCSCData", FALSE)
                ucsc@trackLine <- line
              }
              ucsc
            }
            tsets <- lapply(seq_along(trackLines), makeTrackSet)
            if (drop && length(tsets) == 1)
              tsets[[1]]
            else if (asRangedData)
              do.call(RangedDataList, tsets)
            else
              do.call(GRangesList, tsets)
          })

############ INTERNAL API ############

## every cgi variable is stored in the 'cart'
setClass("ucscCart", contains = "character")

setGeneric("ucscCart", function(object, ...) standardGeneric("ucscCart"))

setMethod("ucscCart", "UCSCSession",
          function(object, form = ucscForm(activeView(object)))
          {
            node <- ucscGet(object, "cart", form)
            raw <- xmlValue(getNodeSet(node, "//pre/text()")[[1]])
            lines <- strsplit(raw, "\n")[[1]]
            fields <- strsplit(lines, " ")
            pairs <- fields[sapply(fields, length) == 2]
            mat <- matrix(unlist(pairs), nrow = 2)
            vals <- mat[2,]
            names(vals) <- mat[1,]
            new("ucscCart", vals)
          })
setMethod("ucscCart", "UCSCView",
          function(object)
          {
            ucscCart(object@session, ucscForm(object))
          })

setMethod("range", "ucscCart",
          function(x, ..., na.rm)
          {
            pos <- x["position"]
            posSplit <- strsplit(pos, ":")[[1]]
            range <- as.numeric(gsub(",", "", strsplit(posSplit[2], "-")[[1]]))
            GRangesForUCSCGenome(x[["db"]], posSplit[1],
                                 IRanges(range[1], range[2]))
          })

### track information

setClass("ucscTracks",
         representation(ids = "character", modes = "character"))

setGeneric("ucscTracks", function(object, ...) standardGeneric("ucscTracks"))

setMethod("ucscTracks", "UCSCSession",
          function(object, form = list())
          {
            tracks <- ucscGet(object, "tracks", form)
            nodes <- getNodeSet(tracks, "//select/option[@selected]/text()")
            trackModes <- sapply(nodes, xmlValue)
            nodes <- getNodeSet(tracks, "//select/@name")
            trackIds <- unlist(nodes)
            ##trackIds <- sapply(nodes, xmlValue)
            nodes <- getNodeSet(tracks, "//select/../a/text()")
            names(trackIds) <- sub("^ ", "", sapply(nodes, xmlValue))
            new("ucscTracks", ids = trackIds, modes = trackModes)
          })

setMethod("ucscTracks", "UCSCView",
          function(object)
          {
            ucscTracks(object@session, ucscForm(object))
          })

setMethod("ucscTrackModes", "ucscTracks",
          function(object)
          {
            modes <- object@modes
            labels <- names(object@ids)
            names(modes) <- object@ids
            ucscTrackModes(modes, labels)
          })

## List available UCSC genomes

ucscGenomes <- function() {
  doc <- httpGet("http://genome.ucsc.edu/FAQ/FAQreleases")
  table <- getNodeSet(doc, "//table[@border='1']")[[1]]
  species <- sapply(getNodeSet(table, "tr/td[1]//text()"), xmlValue)
  species <- sub("^ *", "", species) # attempt to strip weird characters
  # The code below tries to detect the empty cells in the SPECIES col of
  # the table. The real content of these cells seems to vary from one
  # platform to the other (not clear why, maybe some sort of local issue?).
  # There must be a simplest way.
  # TODO: Test this on Windows!
  is_empty_species <- species %in% c("<c2><a0>", "\xc2\xa0", "\xc3\x82\xc2\xa0")
  species <- rep.int(species[!is_empty_species],
                     diff(which(c(!is_empty_species, TRUE))))
  dbs <- sapply(getNodeSet(table, "tr/td[2]//text()"), xmlValue)
  dates <- sapply(getNodeSet(table, "tr/td[3]//text()"), xmlValue)
  nms <- sapply(getNodeSet(table, "tr/td[4]//text()"), xmlValue)
  df <- data.frame(db = dbs, organism = species, date = dates, name = nms)
  status <- getNodeSet(table, "tr/td[5]//text()")
  df <- df[sapply(status, xmlValue) == "Available",]
  rownames(df) <- NULL
  df
}

# form creation

setGeneric("ucscForm", function(object, ...) standardGeneric("ucscForm"))

setMethod("ucscForm", "RangesList",
          function(object)
          {
            form <- list()
            if (length(genome(object)))
              form <- c(form, db = genome(object))
            chrom <- chrom(object)
            if (!is.null(chrom)) {
              if (!length(chrom))
                chrom <- levels(chrom)[1]
              scipen <- getOption("scipen")
              options(scipen = 100) # prevent use of scientific notation
              on.exit(options(scipen = scipen))
              position <- chrom
              if (length(unlist(start(object))))
                position <- paste(ucscNormSeqNames(chrom), ":",
                                  unlist(start(object)), "-",
                                  unlist(end(object)), sep = "")
              form <- c(form, position = position)
            }
            form
          })
setMethod("ucscForm", "GRanges",
          function(object)
          {
            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))
            form <- list()
            if (length(genome(object)))
              form <- c(form, db = genome(object))
            object <- object[1]
            c(form, position = paste(ucscNormSeqNames(seqnames(object)), ":",
                      unlist(start(object)), "-",
                      unlist(end(object)), sep = ""))
          })

setMethod("ucscForm", "UCSCTrackModes",
          function(object)
          {
            as.list(object)
          })
setMethod("ucscForm", "UCSCView",
          function(object)
          {
            if (length(object@hgsid))
              list(hgsid = as.character(object@hgsid))
            else list()
          })
setMethod("ucscForm", "RangedDataList",
          function(object, format, ...)
          {
            lines <- export(object, format = "ucsc", subformat = format, ...)
            text <- paste(paste(lines, collapse = "\n"), "\n", sep = "")
            filename <- paste("track", format, sep = ".")
            upload <- fileUpload(filename, text, "text/plain")
            form <- list(Submit = "Submit", hgt.customFile = upload)
            genome <- genome(object[[1]])
            if (length(genome))
              form <- c(form, db = genome)
            form
          })

setMethod("ucscForm", "UCSCTableQuery",
          function(object, tracks = FALSE) {
            ## range (ie genome) is required
            range <- object@range
            form <- ucscForm(range)
            if (is.null(object@track) && !tracks)
              form <- c(form, list(hgta_group = "allTables"))
            else
              form <- c(form, list(hgta_group = "allTracks",
                                   hgta_track = object@track))
            if (spansGenome(range))
              regionType <- "genome"
            else regionType <- "range"
            form <- c(form, hgta_regionType = regionType)
            table <- object@table
            form <- c(form, hgta_table = table)
            if (!is.null(object@outputType)) {
              form <- c(form, hgta_outputType = object@outputType)
            }
            form
          })

setMethod("ucscForm", "NULL", function(object) list())

# Transforming to a cookie

cookiePair <- function(key, value) paste(key, value, sep = "=")

setGeneric("ucscCookie", function(object, ...) standardGeneric("ucscCookie"))
setMethod("ucscCookie", "UCSCSession",
          function(object)
          {
            cookiePair("hguid", object@hguid)
          })

# HTTP wrappers

# URL constants for UCSC
ucscURLTable <- c(gateway = "hgGateway", tracks = "hgTracks",
                  custom = "hgCustom", tables = "hgTables",
                  cart = "cartDump")

ucscURL <-
  function(object, key)
  {
    path <- ucscURLTable[key]
    if (is.na(path))
      stop("Key '", key, "' does not match a known URL")
    else paste(object@url, path, sep="")
  }

# convenience wrappers for _initialized_ sessions
ucscShow <- function(object, key, .form = list(), ...)
  httpShow(ucscURL(object, key), .form, ...)
ucscPost <- function(object, key, .form = list(), ...)
  httpPost(ucscURL(object, key), .form, ..., cookie = ucscCookie(object))
ucscGet <- function(object, key, .form = list(), ...)
  httpGet(ucscURL(object, key), .form, ..., cookie = ucscCookie(object))
