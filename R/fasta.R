### =========================================================================
### FASTA support
### -------------------------------------------------------------------------
###
### Since we have 2bit, we might as well have FASTA
###

setClass("FastaFile", contains = "RTLFile")

FastaFile <- function(resource) {
  new("FastaFile", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setMethod("export", c("ANY", "FastaFile"), function(object, con, ...) {
  export(as(object, "DNAStringSet"), con, ...)
})

setMethod("export", c("BSgenome", "FastaFile"),
          function(object, con, format, ...) {
            append <- FALSE
            for (seqname in seqnames(object)) {
              dna <- object[[seqname]]
              masks(dna) <- NULL
              dna <- DNAStringSet(dna)
              names(dna) <- seqname
              writeXStringSet(dna, path(con), append = append, ...)
              append <- TRUE
            }
          })

setMethod("export", c("XStringSet", "FastaFile"),
          function(object, con, format, ...)
          {
            writeXStringSet(object, path(con), format = "fasta", ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setMethod("import", "FastaFile",
          function(con, format, text, type = c("DNA", "RNA", "AA", "B"), ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            readFun <- get(paste("read.", match.arg(type), "StringSet",
                                 sep = ""))
            readFun(path(con), format = "fasta", ...)
          })
