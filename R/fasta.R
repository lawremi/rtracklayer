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

setGeneric("export.fasta",
           function(object, con, ...) standardGeneric("export.fasta"))

setMethod("export.fasta", "ANY", function(object, con, ...) {
  export(object, con, "fasta", ...)
})

setMethod("export", c("ANY", "FastaFile"), function(object, con, ...) {
  export(as(object, "DNAStringSet"), con, ...)
})

setMethod("export", c("BSgenome", "FastaFile"),
          function(object, con, format) {
            append <- FALSE
            for (seqname in seqnames(object)) {
              dna <- DNAStringSet(object[[seqname]])
              names(dna) <- seqname
              write.XStringSet(dna, path(con), append = append)
              append <- TRUE
            }
          })

setMethod("export", c("DNAStringSet", "FastaFile"),
          function(object, con, format, append = append)
          {
            write.XStringSet(object, path(con), format = "fasta",
                             append = append)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.fasta", function(con, ...) standardGeneric("import.fasta"))

setMethod("import.fasta", "ANY",
          function(con, ...)
          {
            import(con, "fasta", ...)
          })

setMethod("import", "FastaFile",
          function(con, format, text, type = c("DNA", "RNA", "AA", "B"), ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            readFun <- get(paste("read.", match.arg(type), "StringSet",
                                 sep = ""))
            readFun(path(con), format = "fasta", ...)
          })
