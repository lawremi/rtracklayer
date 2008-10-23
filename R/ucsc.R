# UCSC genome browser interface

# every UCSC session is identified by a 'hgsid'
setClass("ucscSession",
         representation(url = "character", hguid = "numeric",
                        views = "environment"),
         contains = "browserSession")

# gets an 'hgsid' to initialize the session
setMethod("initialize", "ucscSession",
          function(.Object, url = "http://genome.ucsc.edu/cgi-bin/", ...)
          {
            .Object@url <- url
            .Object@views <- new.env()
            handle <- getCurlHandle(...)
            getURL(ucscURL(.Object, "gateway"), cookiefile = tempfile(),
                   curl = handle)
            cart <- getURL(ucscURL(.Object, "cart"), curl = handle)
            .Object@hguid <- as.numeric(gsub(".*hguid=([^\n]*).*", "\\1", cart))
            .Object
          })

setMethod("layTrack", c("ucscSession", "trackSets"),
          function(object, track, name = names(track), view,
                   format = c("auto", "bed", "wig", "gff1"), ...) {
            format <- match.arg(format)
            if (length(track@.Data)) {
              ## upload tracks in blocks, one for each genome
              genomes <- lapply(track, slot, "genome")
              genomes[sapply(genomes, length) == 0] <- ""
              names(track) <- name
              tapply(track, unlist(genomes),
                     function(tracks)
                     {
                       form <- ucscForm(trackSets(tracks), format, ...)
                       response <- ucscPost(object, "custom", form)
### FIXME: need to check for error
                     })
              args <- list()
              if (view) { # optionally view the track
                args <- c(args, genomeSegment(tail(track, 1)[[1]]))
                ## update the view
                do.call("browserView", c(object, args))
              }
            }
            object
          })

setMethod("browserViews", "ucscSession",
          function(object) object@views$instances)

# get the list of track names
setMethod("tracks", "ucscSession",
          function(object) ucscTracks(object)@ids)

# get the current genomeSegment
setMethod("genomeSegment", "ucscSession",
          function(object) genomeSegment(ucscCart(object)))

# export data from UCSC (internal utility)
ucscExport <- function(object, segment, track, table, output, followup = NULL)
{
    segment <- merge(genomeSegment(object), segment)
    segForm <- ucscForm(segment)
    if (length(segment@chrom))
      regionType <- "range"
    else regionType <- "genome"
    form <- c(segForm, list(hgta_table = table, hgta_track = track,
                            hgta_outputType = output,
                            hgta_doTopSubmit = "get output",
                            hgta_regionType = regionType))
    output <- ucscGet(object, "tables", form, .parse = !is.null(followup))
    if (!is.null(followup)) {
      node <- getNodeSet(output, "//input[@name = 'hgsid']/@value")[[1]]
      hgsid <- node ##xmlValue(node)
      form <- c(followup, list(hgsid = hgsid))
      output <- ucscGet(object, "tables", form, .parse = FALSE)
    }
    output
}

## download a trackSet by name
setMethod("trackSet", "ucscSession",
          function(object, name, segment = genomeSegment(object), table = NULL)
          {
            trackids <- tracks(object)
            if (!(name %in% trackids)) {
              mapped_name <- trackids[name]
              if (is.na(mapped_name))
                stop("Unknown track: ", name)
              name <- mapped_name
            }
            if (is.null(table))
              table <- name # default table is track id
            followup <- NULL
            tables <- ucscGet(object, "tables")
            types_path <- "//select[@name = 'hgta_outputType']/option/@value"
            ##types <- sapply(getNodeSet(tables, types_path), xmlValue)
            types <- unlist(getNodeSet(tables, types_path))
            if ("wigData" %in% types) { # track stored as wig
              format <- "wig"
              output <- "wigData"
            } else {
              format <- output <- "bed"
              followup <- list(hgta_doGetBed = "get BED",
                               hgta_printCustomTrackHeaders = "on",
                               boolshad.hgta_printCustomTrackHeaders = "1")
            }
            output <- ucscExport(object, segment, name, table, output, followup)
            import(text = output, format = format)
          })

## grab sequences for features in 'track' at 'segment'
setMethod("genomeSequence", "ucscSession",
          function(object, segment, table = "gold")
          {
            followup <- list(hgta_doGenomicDna = "get sequence",
                             hgSeq.casing = "upper",
                             hgSeq.repMasking = "lower")
            output <- ucscExport(object, segment, "gold", table, "sequence",
                                 followup)
            con <- file()
            writeLines(output, con)
            set <- read.DNAStringSet(con, "fasta")
            close(con)
            set
          })

## get a data.frame from a UCSC table
setGeneric("ucscTable",
           function(object, segment, track, table) standardGeneric("ucscTable"))
setMethod("ucscTable", "ucscSession",
          function(object, segment, track, table)
          {
            output <- ucscExport(object, segment, track, table, "primaryTable")
            read.table(output, sep = "\t")
          })

## UCSC genome view
setClass("ucscView", representation(hgsid = "numeric"),
         contains = "browserView")

## create a view for the given session, position and track visibility settings
## if 'tracks' is a character vector (but not a ucscTrackModes instance) it is
## assumed to name the tracks that should be in the view. otherwise, an
## attempt is made to coerce it to a ucscTrackModes instance.
setMethod("browserView", "ucscSession",
          function(object, segment, track, ...)
          {
            view <- new("ucscView", session = object)
            form <- list()
            seg <- NULL                            # figure out segment
            args <- list(...)
            argsForSeg <- names(args) %in% slotNames("genomeSegment")
            segArgs <- args[argsForSeg]
            if (length(segArgs)) {
              if (missing(segment))
                seg <- do.call("genomeSegment", segArgs)
              else seg <- do.call("genomeSegment", c(segArgs, segment=segment))
            } else if (!missing(segment))
              seg <- segment
            if (!is.null(seg)) # only need to pass it if specified
              form <- c(form, ucscForm(segment))
            ## figure out track modes
            modes <- NULL
            if (!missing(track)) {
              if (is(track, "ucscTrackModes"))
                modes <- track
              else if (class(track) == "character") {
                modes <- ucscTrackModes(object)
                tracks(modes) <- track
                modes <- modes[track]
              } else modes <- as(track, "ucscTrackModes")
            }
            ## new hgsid for each browser launch
            doc <- ucscGet(object, "gateway")
            node <- getNodeSet(doc, "//input[@name = 'hgsid']/@value")[[1]]
            hgsid <- node ##xmlValue(node)
            view@hgsid <- as.numeric(hgsid)
            argModes <- do.call("ucscTrackModes", args[!argsForSeg])
            if (is.null(modes)) # obviously inefficient through here...
              modes <- ucscTrackModes(view)[names(argModes)]
            modes[names(argModes)] <- argModes
            form <- c(form, ucscForm(modes), ucscForm(view))
            ## launch a web browser
            ucscShow(object, "tracks", form)
            ## save this view
            object@views$instances <- c(object@views$instances, view)
            view
          })

# every view has a "mode" (hide, dense, pack, squish, full) for each track
### FIXME: probably should be merged with ucscTracks
### FIXME: and maybe hide the structure entirely, using [ on ucscView
setClass("ucscTrackModes", representation(labels = "character"),
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
            new("ucscTrackModes", object, labels = as.character(labels))
          })
setMethod("ucscTrackModes", "missing",
          function(object, ...) ucscTrackModes(character(), ...))

### FIXME: the cart is not reliable, need to parse hgTracks page

setMethod("ucscTrackModes", "ucscView",
          function(object)
          {
            ucscTrackModes(ucscTracks(object))
          })

setMethod("ucscTrackModes", "ucscSession",
          function(object)
          {
            ucscTrackModes(ucscTracks(object))
          })

setGeneric("ucscTrackModes<-",
           function(object, value) standardGeneric("ucscTrackModes<-"))
setReplaceMethod("ucscTrackModes", c("ucscView", "ucscTrackModes"),
                 function(object, value)
                 { # FIXME: needs to be more extensible
                   browserView(object@session, genomeSegment(object), value)
                 })
setReplaceMethod("ucscTrackModes", c("ucscView", "character"),
                 function(object, value)
                 {
                   ucscTrackModes(object) <- ucscTrackModes(value)
                   object
                 })

## subsetting ucscTrackModes

## if not in ids, try labels
resolveTrackIndex <- function(object, i) {
  if (is.character(i)) {
    missing <- !(i %in% names(object))
    matching <- match(i[missing], object@labels)
    if (any(is.na(matching))) {
      unmatched <- i[missing][is.na(matching)]
      stop("Unknown tracks:", paste(unmatched, collapse = ", "))
    }
    i[missing] <- names(object)[matching]
  }
  i
}

setMethod("[", "ucscTrackModes", function(x, i, j, ..., drop=FALSE) {
  vec <- x@.Data
  names(vec) <- names(x)
  names(x@labels) <- names(x)
  ind <- resolveTrackIndex(x, i)
  initialize(x, vec[ind], labels = x@labels[ind])
})

setReplaceMethod("[", "ucscTrackModes", function(x, i, j, ..., value) {
  vec <- x@.Data
  names(vec) <- names(x)
  vec[resolveTrackIndex(x, i)] <- value
  x@.Data <- as.vector(vec)
  x
})

# handle simple track show/hide

setMethod("tracks", "ucscTrackModes",
          function(object)
          {
            visible <- object != "hide"
            tracks <- names(object)[visible]
            names(tracks) <- object@labels[visible]
            tracks
          })
setReplaceMethod("tracks", "ucscTrackModes",
                 function(object, value)
                 {
                   value <- resolveTrackIndex(object, value)
                   spec <- names(object) %in% value
                   object[!spec] <- "hide"
                   object[spec & object == "hide"] <- "full"
                   object
                 })

setMethod("tracks", "ucscView",
          function(object)
          {
            tracks <- ucscTracks(object)
            modes <- ucscTrackModes(tracks)
            tracks@ids[tracks@ids %in% tracks(modes)]
          })
setReplaceMethod("tracks", "ucscView",
                 function(object, value)
                 {
                   tracks(ucscTrackModes(object)) <- value
                   object
                 })

setMethod("genomeSegment", "ucscView",
          function(object) genomeSegment(ucscCart(object)))
setReplaceMethod("genomeSegment", "ucscView",
                 function(object, value)
                 {
                   # need to check for partially specified segment
                   # and resolve that here
                   browserView(object@session, value, ucscTrackModes(object))
                 })

# only one view per session; a view is always active
setMethod("activeView", "ucscView", function(object) TRUE)

# ucscTrackSet

# visual properties are specified by a "track line" for UCSC
setClass("ucscTrackLine",
         representation(name = "character", description = "character",
                        visibility = "character", color = "integer",
                        priority = "numeric"),
         prototype(name = "R Track"))

setMethod("show", "ucscTrackLine",
          function(object)
          {
            cat(as(object, "character"), "\n")
          })

setClass("basicTrackLine",
         representation(itemRgb = "logical", useScore = "logical",
                        group = "character", db = "character",
                        offset = "numeric", url = "character",
                        htmlUrl = "character"),
         contains = "ucscTrackLine")

ucscPair <- function(key, value) paste(key, value, sep = "=")

# to a text line
setAs("ucscTrackLine", "character",
      function(from)
      {
        checkString <- function(str, len) {
          if (nchar(gsub("[a-zA-Z0-9_ ]", "", str)))
            warning("The string '", str,
                    "' contains non-standard characters.")
          if (nchar(str) > len) {
            str <- substring(str, 1, len)
            warning("The string '", str, "' must be less than ", len,
                    " characters; it has been truncated.")
          }
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

setAs("character", "basicTrackLine",
      function(from)
      {
        str <- as(as(from, "ucscTrackLine"), "character")
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
setAs("character", "ucscTrackLine",
      function(from)
      {
        line <- new("ucscTrackLine")
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

setAs("character", "basicTrackLine",
      function(from)
      {
        line <- new("basicTrackLine", as(from, "ucscTrackLine"))
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

# each track in UCSC has a "track line"
# other software also understands this (i.e. IGB)
setClass("ucscTrackSet",
         representation(trackLine = "ucscTrackLine"),
         prototype(trackLine = new("basicTrackLine")),
         "trackSet")

setMethod("initialize", "ucscTrackSet",
          function(.Object,
                   assayData = assayDataNew(
                     dataVals = as.matrix(dataVals), ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   dataVals = matrix(nrow = ifelse(missing(featureData), 0,
                                       nrow(pData(featureData))), ncol = 0),
                   genome = "hg18",
                   trackLine = new("basicTrackLine"),
                   ...)
          {
            if (!missing(assayData) && !missing(dataVals))
              stop("'assayData' and 'dataVals' cannot both be specified")
            .Object@trackLine <- trackLine
            callNextMethod(.Object,
                           assayData = assayData,
                           phenoData = phenoData,
                           featureData = featureData,
                           experimentData = experimentData,
                           annotation = annotation,
                           genome = genome)
          })

setMethod("show", "ucscTrackSet",
          function(object)
          {
            callNextMethod()
            cat("trackLine:", as(object@trackLine, "character"), "\n")
          })

# the 'ucsc' format is a meta format with a track line followed by
# tracks formatted as 'wig', 'bed', 'gff', 'gtf', or 'psl'.
# currently, only gff and wig are supported
setGeneric("export.ucsc",
           function(object, con, subformat = c("auto", "gff1", "wig", "bed"),
                    ...)
           standardGeneric("export.ucsc"))

setMethod("export.ucsc", "trackSets",
          function(object, con, subformat = c("auto", "gff1", "wig", "bed"),
                   trackNames, ...)
          {
            subformat <- match.arg(subformat)
            if (missing(trackNames)) {
              trackNames <- names(object)
              if (is.null(trackNames))
                trackNames <- paste("R Track", seq_along(object))
              ucsc <- sapply(object, is, "ucscTrackSet")
              lines <- sapply(object[ucsc], slot, "trackLine")
              trackNames[ucsc] <- as.character(sapply(lines, slot, "name"))
            }
            for (i in seq_along(object))
              export.ucsc(object[[i]], con, subformat, name = trackNames[i],
                          ...)
          })

# wig is a special case (requires a special track line)
ucscTrackLineClass <- function(subformat)
{
  if (subformat == "wig")
    "wigTrackLine"
  else "basicTrackLine"
}

setMethod("export.ucsc", "trackSet",
          function(object, con, subformat, ...)
          {
            object <- as(object, "ucscTrackSet")
            subformat <- match.arg(subformat)
            export.ucsc(object, con, subformat, ...)
          })
setMethod("export.ucsc", "ucscTrackSet",
          function(object, con, subformat, ...)
          {
            subformat <- match.arg(subformat)
            auto <- FALSE
            if (subformat == "auto") {
              auto <- TRUE
              subformat <- "bed"
              if (is.numeric(dataVals(object)))
                subformat <- "wig"
            }
            if (subformat == "wig") {
              strand <- as.character(strand(object))
              strand[is.na(strand)] <- "NA"
              isDisjoint <- function(track) {
                starts <- start(track)
                ends <- end(track)
                startord <- order(starts)
                all(tail(starts[startord], -1) - head(ends[startord], -1) > 0)
              }
              if (!all(sapply(split(object, strand), isDisjoint))) {
                if (auto)
                  subformat <- "bed"
                else stop("Track not compatible with WIG format: ",
                          "Overlapping features must be on separate strands")
              }
            }
            lineClass <- ucscTrackLineClass(subformat)
            if (!is(object@trackLine, lineClass))
              object@trackLine <- as(object@trackLine, lineClass)
            args <- list(...)
            lineArgs <- names(args) %in% slotNames(lineClass)
            for (argName in names(args)[lineArgs])
              slot(object@trackLine, argName) <- args[[argName]]
            if (subformat == "wig") {
              subformat <- "wigLines"
              strand <- as.character(strand(object))
              strand[is.na(strand)] <- "NA"
              if (!all(strand[1] == strand)) {
                nameMap <- c("+" = "p", "-" = "m", "NA" = "NA")
                strand <- factor(strand)
                name <- paste(object@trackLine@name, nameMap[levels(strand)])
                tracks <- split(object, strand)
                export.ucsc(tracks, con, "wig", name, ...)
                return()
              }
            }
            writeLines(as(object@trackLine, "character"), con)
            do.call("export", c(list(as(object, "trackSet"), con, subformat),
                                args[!lineArgs]))
          })

# for GFF, the track line should go in a comment
setMethod("export.gff", "ucscTrackSet",
          function(object, con, version, source)
          {
            gffComment(con, as(object@trackLine, "character"))
            callNextMethod()
          })

setGeneric("import.ucsc",
           function(con, subformat = c("auto", "gff1", "wig", "bed"),
                    drop = FALSE, ...)
           standardGeneric("import.ucsc"))
setMethod("import.ucsc", "ANY",
          function(con, subformat, drop = FALSE, ...)
          {
            subformat <- match.arg(subformat)
            lines <- readLines(con, warn = FALSE)
            tracks <- grep("^track", lines)
            trackLines <- lines[tracks]
            starts <- tracks+1
            ends <- c(tail(tracks,-1)-1, length(lines))
            makeTrackSet <- function(i)
            {
              line <- as(trackLines[i], ucscTrackLineClass(subformat))
              if (subformat == "wig")
                subformat <- "wigLines"
              text <- lines[starts[i]:ends[i]]
              trackSet <- import(format = subformat, text = text, ...)
              ucsc <- as(trackSet, "ucscTrackSet", FALSE)
              ucsc@trackLine <- line
              ucsc
            }
            tsets <- lapply(seq_along(trackLines), makeTrackSet)
            if (drop && length(tsets) == 1)
              tsets[[1]]
            else trackSets(tsets)
          })

############ INTERNAL API ############

# every cgi variable is stored in the 'cart'
setClass("ucscCart", contains = "character")

setGeneric("ucscCart", function(object, ...) standardGeneric("ucscCart"))

setMethod("ucscCart", "ucscSession",
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
setMethod("ucscCart", "ucscView",
          function(object)
          {
            ucscCart(object@session, ucscForm(object))
          })

setMethod("genomeSegment", "ucscCart",
          function(object)
          {
            pos <- object["position"]
            posSplit <- strsplit(pos, ":")[[1]]
            range <- as.numeric(strsplit(posSplit[2], "-")[[1]])
            genomeSegment(genome = object["db"], chrom = posSplit[1],
                          start = range[1], end = range[2])
          })
            
#setMethod("ucscTrackModes", "ucscCart",
#          function(object)
#          {
#            modes <- character()
#            for (mode in c("hide", "dense", "pack", "squish", "full"))
#              modes[names(object)[object == mode]] <- mode
#            ucscTrackModes(modes)
#          })


# track information

setClass("ucscTracks",
         representation(ids = "character", modes = "character"))

setGeneric("ucscTracks", function(object, ...) standardGeneric("ucscTracks"))

setMethod("ucscTracks", "ucscSession",
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

setMethod("ucscTracks", "ucscView",
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

# form creation

setGeneric("ucscForm", function(object, ...) standardGeneric("ucscForm"))

setMethod("ucscForm", "genomeSegment",
          function(object)
          {
            form <- list()
            if (length(object@genome))
              form <- c(form, db = object@genome)
            if (length(object@chrom)) {
              scipen <- getOption("scipen")
              options(scipen = 100) # prevent use of scientific notation
              on.exit(options(scipen = scipen))
              position <- paste(object@chrom, ":", object@start, "-",
                                object@end, sep = "")
              form <- c(form, position = position)
            }
            form
          })
setMethod("ucscForm", "ucscTrackModes",
          function(object)
          {
            as.list(object)
          })
setMethod("ucscForm", "ucscView",
          function(object)
          {
            if (length(object@hgsid))
              list(hgsid = as.character(object@hgsid))
            else list()
          })
setMethod("ucscForm", "trackSets",
          function(object, format, ...)
          {  
            lines <- export(object, format = "ucsc", subformat = format, ...)
            text <- paste(paste(lines, collapse = "\n"), "\n", sep = "")
            filename <- paste("track", format, sep = ".")
            upload <- fileUpload(filename, text, "text/plain")
            form <- list(Submit = "Submit", hgt.customFile = upload)
            genome <- object[[1]]@genome
            if (length(genome))
              form <- c(form, db = genome)
            form
          })

setMethod("ucscForm", "NULL", function(object) list())

# Transforming to a cookie

cookiePair <- function(key, value) paste(key, value, sep = "=")

setGeneric("ucscCookie", function(object, ...) standardGeneric("ucscCookie"))
setMethod("ucscCookie", "ucscSession",
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
