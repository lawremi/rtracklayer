# central abstraction, represents a genome browser
setClass("BrowserSession")

# a single genome view in a session
setClass("BrowserView", representation(session = "BrowserSession"),
         contains = "VIRTUAL")

# create a browser view
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

# get/set (names of) tracks from e.g. a view
setGeneric("trackNames", function(object, ...) standardGeneric("trackNames"))
setGeneric("trackNames<-",
           function(object, value) standardGeneric("trackNames<-"))

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
                paste(names(range), ":", start(range), "-", end(range),
                      sep = ""),
                "\n")
            cat(IRanges:::labeledLine("trackNames", names(trackNames(object))))
          })

setGeneric("track<-",
           function(object, ..., value) standardGeneric("track<-"))
# load a track into a browser
setReplaceMethod("track", c("BrowserSession", "RangedData"),
          function(object, name = deparse(substitute(value)), view = FALSE, ...,
                   value)
                 {
                   track(object, name, view, ...) <- RangedDataList(value)
                   object
                 })
# load several tracks into a browser
# (this may be more efficient for some implementations)
setReplaceMethod("track", c("BrowserSession", "RangedDataList"),
                 function(object, name = names(track), view = FALSE, ..., value)
                 {
                   for (i in seq_len(length(name) - 1))
                     track(object, name[i], FALSE, ...) <- value[[i]]
                   last <- tail(value, 1)
                   if (length(last))
                     track(object, tail(name, 1), view, ...) <- last[[1]]
                   object
                 })

setClassUnion("RangedDataORRangedDataList", c("RangedData", "RangedDataList"))

setMethod("[[<-", c("BrowserSession", value="RangedDataORRangedDataList"),
          function(x, i, j, ..., value) {
            if (!missing(j))
              warning("argument 'j' ignored")
            track(x, i, ...) <- value
            x
          })

setMethod("$<-", c("BrowserSession", value="RangedDataORRangedDataList"),
          function(x, name, value) {
            x[[name]] <- value
            x
          })

setGeneric("track", function(object, ...) standardGeneric("track"))

setMethod("[[", "BrowserSession", function (x, i, j, ...) {
  if (!missing(j))
    warning("argument 'j' ignored")
  track(x, i, ...)
})

setMethod("$", "BrowserSession", function (x, name) {
  x[[name]]
})

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

# high-level entry point
setGeneric("browseGenome",
           function(object, ...)
           standardGeneric("browseGenome"))

setMethod("browseGenome", "missing",
          function(object, ...) browseGenome(RangedDataList(), ...))

setMethod("browseGenome", "RangedDataORRangedDataList",
          function(object, browser = "UCSC",
                   range = base::range(object), view = TRUE,
                   trackParams = list(), viewParams = list(), ...)
          {
            # initialize session of type identified by 'browser'
            session <- browserSession(browser, ...)
            # load 'object'
            trackParams <- c(list(session), trackParams)
            pcall <- sys.call(sys.parent(1))
            name <- deparse(as.list(match.call(call=pcall))[[2]])
            if (is(object, "RangedData"))
              trackParams <- c(trackParams, name = name)
            session <- do.call(`track<-`, c(trackParams, list(value = object)))
            # open view of 'range'
            if (view) {
              range <- mergeRange(range(session), range)
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

# load a sequence into the browser
setGeneric("sequence<-", function(object, name, ..., value)
           standardGeneric("sequence<-"))
