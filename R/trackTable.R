### =========================================================================
### Arbitrary track table support
### -------------------------------------------------------------------------
###
### At this point, this is internal, in support of the tabix stuff
###

setGeneric("import.trackTable",
           function(con, genome = NA, seqnames = 1L, start = 2L, end = 3L, ...)
           standardGeneric("import.trackTable"),
           signature = "con")

setMethod("import.trackTable", "characterORconnection",
          function(con, genome = NA, seqnames = 1L, start = 2L, end = 3L, ...)
          {
            tab <- read.table(con, sep = "\t", ...)
            GenomicData(IRanges(tab[[start]], tab[[end]]),
                        tab[-c(seqnames, start, end)],
                        chrom = tab[[seqnames]], genome = genome)
          })

