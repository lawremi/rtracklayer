### =========================================================================
### Arbitrary track table support
### -------------------------------------------------------------------------
###
### At this point, this is internal, in support of the tabix stuff
###

setClass("TabSeparatedFile", contains = "BiocFile")
TabSeparatedFile <- function(resource) {
  new("TabSeparatedFile", resource = resource)
}

setGeneric("import.tabSeparated",
           function(con, genome = NA, seqnames = 1L, start = 2L, end = 3L, ...)
           standardGeneric("import.tabSeparated"),
           signature = "con")

setMethod("import.tabSeparated", "character_OR_connection",
          function(con, genome = NA, seqnames = 1L, start = 2L, end = 3L, ...)
          {
            tab <- read.table(con, sep = "\t", ...)
            ans <- GRanges(tab[[seqnames]], IRanges(tab[[start]], tab[[end]]),
                           strand=Rle(strand("*"), nrow(tab)),
                           tab[-c(seqnames, start, end)])
            metadata(ans) <- list(genome = genome)
            ans
          })

setGeneric("export.tabSeparated",
           function(object, con, ...) standardGeneric("export.tabSeparated"))

setMethod("export.tabSeparated", "ANY",
          function(object, con, ...) {
            export(object, con, "tabSeparated", ...)
          })

setMethod("export", c("ANY", "TabSeparatedFile"),
          function(object, con, ...) {
            df <- as.data.frame(object)
            write.table(df, path(con), sep="\t", ...)
          })

setMethod("export", c("GenomicRanges", "TabSeparatedFile"),
          function(object, con, ...) {
            object <- as.data.frame(object, row.names=NULL)
            object$width <- NULL
            export(object, con, ...)
          })
