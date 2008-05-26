# central abstraction, represents a genome browser
setClass("browserSession")

# a single genome view in a session
setClass("browserView", representation(session = "browserSession"),
         contains = "VIRTUAL")

# create a browser view
setGeneric("browserView",
           function(object, segment = genomeSegment(object),
                    track = tracks(object), ...)
           standardGeneric("browserView"))

# get the views from the browser
setGeneric("browserViews",
           function(object, ...) standardGeneric("browserViews"))

# get/set active view
setGeneric("activeView", function(object) standardGeneric("activeView"))
setGeneric("activeView<-", function(object, value)
           standardGeneric("activeView<-"))

setMethod("activeView", "browserSession",
          function(object)
          {
            views <- browserViews(object)
            active <- sapply(views, activeView)
            if (any(active))
              views[[tail(which(active),1)]]
            else NULL
          })

setMethod("show", "browserSession",
          function(object)
          {
            cat("A genome browser session of class '", class(object), "' with ",
                length(browserViews(object)), " views and ",
                length(tracks(object)), " tracks\n", sep ="")
          })

# close a session or view
setGeneric("close", function(con, ...) standardGeneric("close"))

# FIXME: what about isOpen?

# get/set (names of) tracks from e.g. a view
setGeneric("tracks", function(object, ...) standardGeneric("tracks"))
setGeneric("tracks<-", function(object, value) standardGeneric("tracks<-"))

# get/set the selected features in e.g. a view
# this can return a list, with a logical vector for each track
setGeneric("selectedFeatures", function(object, ...)
           standardGeneric("selectedFeatures"))
setGeneric("selectedFeatures<-", function(object, value)
           standardGeneric("selectedFeatures<-"))

setMethod("show", "browserView", function(object)
          {
            cat("A genome browser view of class '", class(object), "'.\n\n",
                sep = "")
            cat("Segment:\n")
            show(genomeSegment(object))
            cat("Tracks:\n")
            show(tracks(object))
          })

setGeneric("layTrack",
           function(object, track, name = deparse(substitute(track)),
                    view = FALSE, ...)
           standardGeneric("layTrack"))
# load a track into a browser
setMethod("layTrack", c("browserSession", "trackSet"),
          function(object, track, name, view, ...)
          {
            layTrack(object, trackSets(track), name, view, ...)
          })
# load several tracks into a browser
# (this may be more efficient for some implementations)
setMethod("layTrack", c("browserSession", "trackSets"),
          function(object, track, name = names(track), view, ...)
          {
            for (i in seq_len(length(name) - 1))
            #for (t in head(track, -1))
              layTrack(object, name[i], FALSE, ...) <- track[[i]]
            last <- tail(track, 1)
            if (length(last))
              object <- layTrack(object, last[[1]], tail(name, 1), view, ...)
            object
          })

# get genome segment of active view (or default if no views)
setMethod("genomeSegment", "browserSession",
          function(object, ...)
          {
            view <- activeView(object)
            if (!is.null(view))
              genomeSegment(view, ...)
            else NULL
          })

# high-level entry point
setGeneric("browseGenome",
           function(tracks = trackSets(), browser = "ucsc",
                    segment = genomeSegment(tracks), view = TRUE,
                    trackParams = list(), viewParams = list(), ...)
           standardGeneric("browseGenome"))

setMethod("browseGenome", "ANY",
          function(tracks, browser, segment, view, trackParams, viewParams, ...)
          {
            # initialize session of type identified by 'browser'
            session <- browserSession(browser)
            # load 'tracks'
            trackParams <- c(list(session, tracks, view = FALSE), trackParams)
            session <- do.call("layTrack", trackParams)
            # open view of 'segment'
            if (view) {
              segment <- merge(genomeSegment(session), segment)
              if (length(list(...)))
                segment <- genomeSegment(..., segment = segment)
              viewParams <- c(list(session, segment), viewParams)
              do.call("browserView", viewParams)
            }
            session
          })

# list names of available genome browsers
genomeBrowsers <- function(where = topenv(parent.frame()))
{
  cl <- getClass("browserSession", where = where)
  browsers <- names(cl@subclasses)
  browsers <- browsers[!sapply(browsers, isVirtualClass, where)]
  sub("Session$", "", browsers)
}

# obtain a browser session
setGeneric("browserSession",
           function(object, ...) standardGeneric("browserSession"))

# instantiate a browser session by name
setMethod("browserSession", "character",
          function(object, ...)
          {
            class <- paste(object, "Session", sep = "")
            if (!extends(class, "browserSession"))
              stop("Browser named '", object, "' is unsupported.")
            new(class, ...)
          })

setMethod("browserSession", "missing",
          function(object, ...) browserSession("ucsc", ...))

# get one from a view
setMethod("browserSession", "browserView", function(object) object@session)

# load a sequence into the browser
setGeneric("laySequence", function(object, sequence, name, ...)
           standardGeneric("laySequence"))

# retrieve a segment of a genome sequence from a browser
setGeneric("genomeSequence",
           function(object, segment = genomeSegment(object), ...)
           standardGeneric("genomeSequence"))
