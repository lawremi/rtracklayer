### =========================================================================
### TrackDb class
### -------------------------------------------------------------------------
###
### Represents a database or source of tracks, i.e, interval datasets
###

setClass("TrackDb")

setGeneric("track<-",
           function(object, ..., value) standardGeneric("track<-"))
## load a track into a database
setReplaceMethod("track", c("TrackDb", "ANY"),
                 function(object, name = deparse(substitute(value)), ..., value)
                 {
                   track(object, name, ...) <-
                     RangedDataList(as(value, "RangedData", strict = FALSE))
                   object
                 })

setClassUnion("RangedDataORRangedDataList", c("RangedData", "RangedDataList"))

setMethod("[[<-", c("TrackDb", value="ANY"),
          function(x, i, j, ..., value) {
            if (!missing(j))
              warning("argument 'j' ignored")
            track(x, i, ...) <- value
            x
          })

setMethod("$<-", c("TrackDb", value="ANY"),
          function(x, name, value) {
            x[[name]] <- value
            x
          })

setGeneric("track", function(object, ...) standardGeneric("track"))

setMethod("[[", "TrackDb", function (x, i, j, ...) {
  if (!missing(j))
    warning("argument 'j' ignored")
  track(x, i, ...)
})

setMethod("$", "TrackDb", function (x, name) {
  x[[name]]
})
