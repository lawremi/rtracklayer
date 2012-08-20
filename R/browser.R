## central abstraction, represents a genome browser

## Question: is a browser session a database of tracks? No question
## that a browser has a track database. The question is of
## inheritance. Historically, browser session has had all the behavior
## we expect from a track database. But it is so much more. In
## particular, it has a list of views. We want to separate the data
## from the views. It is tempting to have the model (tracks) and views
## in separate objects, composed together in a session. However, the
## reality is that browser sessions are monolithic: UCSC can view only
## data UCSC knows, IGB can view only data IGB knows. Thus, we use
## simple inheritance here, saving many refactoring headaches.

setClass("BrowserSession", contains = c("TrackDb", "VIRTUAL"))

## alias for names() for clarity in the session context
setGeneric("trackNames", function(object, ...) standardGeneric("trackNames"))
setGeneric("trackNames<-", function(object, ..., value)
           standardGeneric("trackNames<-"))

setMethod("names", "BrowserSession", function(x) trackNames(x))
setMethod("trackNames", "BrowserSession", function(object) names(object))

## a single genome view in a session
setClass("BrowserView", representation(session = "BrowserSession"),
         contains = "VIRTUAL")

# create one or more browser views
setGeneric("browserView",
           function(object, range, track, ...)
           standardGeneric("browserView"))

# get the views from the browser
setGeneric("browserViews",
           function(object, ...) standardGeneric("browserViews"))

# get/set active view
setGeneric("activeView", function(object) standardGeneric("activeView"))
setGeneric("activeView<-", function(object, value)
           standardGeneric("activeView<-"))

setMethod("activeView", "BrowserSession",
          function(object)
          {
            views <- browserViews(object)
            active <- sapply(views, activeView)
            if (any(active))
              views[[tail(which(active),1)]]
            else NULL
          })

setMethod("show", "BrowserSession",
          function(object)
          {
            cat("A genome browser session of class '", class(object), "' with ",
                length(browserViews(object)), " views and ",
                length(trackNames(object)), " tracks\n", sep ="")
          })

# close a session or view
setGeneric("close", function(con, ...) standardGeneric("close"))

# FIXME: what about isOpen?

# get/set visibility of view tracks
setGeneric("visible", function(object, ...) standardGeneric("visible"))
setMethod("visible", "BrowserView", function(object) {
  trackNames(browserSession(object)) %in% trackNames(object)
})
          
setGeneric("visible<-",
           function(object, ..., value) standardGeneric("visible<-"))
setReplaceMethod("visible", "BrowserView", function(object, value) {
  trackNames(object) <- names(value)[value]
  object
})

# get/set the selected features in e.g. a view
# this can return a list, with a logical vector for each track
setGeneric("selectedFeatures", function(object, ...)
           standardGeneric("selectedFeatures"))
setGeneric("selectedFeatures<-", function(object, value)
           standardGeneric("selectedFeatures<-"))

setMethod("show", "BrowserView", function(object)
          {
            range <- range(object)
            cat(class(object), "of",
                paste(names(range), ":", unlist(start(range)), "-",
                      unlist(end(range)), sep = ""),
                "\n")
            nms <- paste("'", names(trackNames(object)), "'", sep = "")
            cat(IRanges:::labeledLine("trackNames", nms))
          })

setClass("BrowserViewList", contains = "SimpleList",
         prototype = prototype(elementType = "BrowserView"))

BrowserViewList <- function(...) {
  views <- list(...)
  if (length(views) == 1 && is.list(views[[1L]]))
    views <- views[[1L]]
  if (!all(sapply(views, is, "BrowserView")))
    stop("all elements in '...' must be BrowserView objects")
  IRanges:::newList("BrowserViewList", views)
}

# get genome range of active view (or default if no views)
setMethod("range", "BrowserSession",
          function(x, ..., na.rm)
          {
            if (length(...) > 0)
              stop("arguments in '...' ignored")
            view <- activeView(x)
            if (!is.null(view))
              range(view)
            else NULL
          })

setGeneric("range<-", function(x, ..., value) standardGeneric("range<-"))

## just the genome, for convenience
setMethod("genome", "BrowserSession", function(x) genome(range(x)))
setReplaceMethod("genome", "BrowserSession", function(x, value) {
  if (!isSingleString(value))
    stop("'genome' must be a single string")
  genome(range(x)) <- value
  x
})

# high-level entry point
setGeneric("browseGenome",
           function(object, ...)
           standardGeneric("browseGenome"))

setMethod("browseGenome", "missing",
          function(object, ...) browseGenome(RangedDataList(), ...))

setMethod("browseGenome", "GRanges",
          function(object, ...) browseGenome(as(object, "RangedData"), ...))

setMethod("browseGenome", "RangedDataORRangedDataList",
          function(object, browser = "UCSC",
                   range = base::range(object), view = TRUE,
                   trackParams = list(), viewParams = list(),
                   name = "customTrack", ...)
          {
            # initialize session of type identified by 'browser'
            session <- browserSession(browser, ...)
            # load 'object'
            trackParams <- c(list(session), trackParams)
            if (is(object, "RangedData"))
              trackParams <- c(trackParams, name = name)
            session <- do.call(`track<-`, c(trackParams, list(value = object)))
            # open view of 'range'
            if (view) {
              if (!missing(range))
                range <- normGenomeRange(range, session)
              viewParams <- c(list(session, range), viewParams)
              do.call(browserView, viewParams)
            }
            session
          })

# list names of available genome browsers
genomeBrowsers <- function(where = topenv(parent.frame()))
{
  cl <- getClass("BrowserSession", where = where)
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
            if (!extends(class, "BrowserSession"))
              stop("Browser named '", object, "' is unsupported.")
            new(class, ...)
          })

setMethod("browserSession", "missing",
          function(object, ...) browserSession("UCSC", ...))

# get one from a view
setMethod("browserSession", "BrowserView", function(object) object@session)

# load a sequence into the browser (probably should remove this)
setGeneric("sequence<-", function(object, name, ..., value)
           standardGeneric("sequence<-"))

