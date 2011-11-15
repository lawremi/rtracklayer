### =========================================================================
### FASTA support
### -------------------------------------------------------------------------
###
### Since we have 2bit, we might as well have FASTA
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.fasta",
           function(object, con, ...) standardGeneric("export.fasta"))

setMethod("export.fasta", "ANY", function(object, con, ...) {
  export.fasta(as(object, "DNAStringSet"), con, ...)
})

setMethod("export.fasta", c("BSgenome", "character"),
          function(object, con) {
            append <- FALSE
            for (seqname in seqnames(object)) {
              writeFASTA(list(object[[seqname]]), con, desc = list(seqname),
                         append = append)
              append <- TRUE
            }
            file
          })

setMethod("export.fasta", c("DNAStringSet", "character"), function(object, con)
          {
            write.XStringSet(object, con, format = "fasta")
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###

setGeneric("import.fasta", function(con, ...) standardGeneric("import.fasta"))

setMethod("import.fasta", "connection",
          function(con, ...)
          {
            import.fasta(summary(con)$description, ...)
          })

setMethod("import.fasta", "character",
          function(con, type = c("DNA", "RNA", "AA", "B"), ...)
          {
            readFun <- get(paste("read.", match.arg(type), "StringSet",
                                 sep = ""))
            readFun(con, format = "fasta", ...)
          })
