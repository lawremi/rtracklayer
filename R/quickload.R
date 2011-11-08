### =========================================================================
### Quickload support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quickload class
###

setClass("Quickload", representation(uri = "character"))

uri <- function(x, ...) x@uri

Quickload_contents <- function(x) {
  read.table(contentsFile(x), sep = "\t", col.names = c("dir", "title"),
             colClasses = "character")
}

setMethod("genome", "Quickload", function(x) {
  contents <- Quickload_contents(x)
  as.character(structure(contents$dir, names = contents$title))
})

setMethod("names", "Quickload", genome)

setMethod("length", "Quickload", function(x) length(names(x)))

setMethod("[[", "Quickload", function (x, i, j, ...) {
  if (!missing(j))
    warning("argument 'j' ignored")
  QuickloadGenome(x, i, ...)
})

setMethod("$", "Quickload", function (x, name) {
  QuickloadGenome(x, name)
})

Quickload <- function(uri = "quickload", create = FALSE) {
  if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  if (create) {
    if (uriExists(uri)) {
      message("NOTE: '", uri, "' already exists")
      create <- FALSE
    } ## must create this before calling normURI (requires existence)
    else createResource(uri, dir = TRUE)
  }
  ql <- new("Quickload", uri = normURI(uri))
  if (create)
    createResource(contentsFile(ql))
  ql
}

setAs("character", "Quickload", function(from) Quickload(from))

setMethod("show", "Quickload", function(object) {
  cat(class(object), "repository\nuri:", uri(object), "\n")
  cat(IRanges:::labeledLine("genomes", genome(object)))
})

contentsFile <- function(x) file.path(uri(x), "contents.txt")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### QuickloadGenome class
###

setClass("QuickloadGenome",
         representation(quickload = "Quickload",
                        genome = "character"),
         contains = "TrackDb")

setMethod("uri", "QuickloadGenome",
          function(x) file.path(uri(quickload(x)), genome(x)))

quickload <- function(x, ...) x@quickload

setMethod("genome", "QuickloadGenome", function(x) x@genome)

setMethod("seqinfo", "QuickloadGenome", function(x) {
  genome_info <- read.table(genomeFile(x), sep = "\t",
                            col.names = c("seqnames", "seqlengths"),
                            colClasses = c("character", "integer"))
  Seqinfo(genome_info$seqnames, genome_info$seqlengths,
          genome = rep(genome(x), nrow(genome_info)))
})

setReplaceMethod("seqinfo", "QuickloadGenome", function(x, value)
                 {
                   if (uriIsWritable(genomeFile(x))) {
                     df <- as.data.frame(value)[1]
                     write.table(df, genomeFile(x),
                                 quote = FALSE, col.names = FALSE, sep = "\t")
                   } else stop("Repository is read only; cannot write seqinfo")
                   x
                 })

setMethod("releaseDate", "QuickloadGenome", function(x) {
  sub(".*?_(.*?)_([^_]*)$", "\\1 \\2", genome(x))
})

setMethod("organism", "QuickloadGenome", function(x) {
  gsub("_", " ", sub("(.*?)_.*?_[^_]*$", "\\1", genome(x)))
})

setMethod("length", "QuickloadGenome", function(x) {
  length(names(x))
})

QuickloadGenome_annotFiles <- function(x) {
  xmlChildren(xmlInternalTreeParse(annotsFile(x)))$files
}

setMethod("names", "QuickloadGenome", function(x) {
  emd <- elementMetadata(x)
  structure(sapply(as.character(emd$name), URLdecode),
            names = as.character(emd$title))
})

setMethod("elementMetadata", "QuickloadGenome", function(x) {
  files <- QuickloadGenome_annotFiles(x)
  if (!length(xmlChildren(files)))
    new("DataFrame", nrows = length(x))
  else
    Reduce(function(x, y) merge(as.data.frame(as.list(x)),
                                as.data.frame(as.list(y)), all = TRUE),
           lapply(xmlChildren(files), xmlAttrs))
})

setMethod("toString", "QuickloadGenome", function(x) {
  Quickload_contents(quickload(x))[genome(x),"title"]
})

addGenomeToContents <- function(x, title) {
  contents <- Quickload_contents(quickload(x))
  if (!genome(x) %in% contents$dir) {
    contents <- rbind(contents, data.frame(genome(x), title))
    if (uriIsWritable(contentsFile(quickload(x))))
      write.table(contents, contentsFile(quickload(x)),
                  quote = FALSE, row.names = FALSE, col.names = FALSE,
                  sep = "\t")
    else stop("Repository is read only; cannot add genome to contents")
  } else warning("Genome '", genome(x), "' already in contents; not replaced")
}

setMethod("toString", "BSgenome", function(x) {
  paste(organism(x), provider(x), providerVersion(x))
})

## Not sure where this method should land. I'm sure it will be useful
## for publishing adhoc genomes through Quickload.
setMethod("seqinfo", "DNAStringSet", function(x) {
  x_names <- names(x)
  if (is.null(x_names))
    x_names <- as.character(seq(length(x)))
  Seqinfo(x_names, width(x))
})

QuickloadGenome <- function(quickload, genome, create = FALSE,
                            seqinfo = GenomicRanges::seqinfo(genome),
                            title = toString(genome))
{
  if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  quickload <- as(quickload, "Quickload")
  genome <- normArgGenome(genome)
  genome_id <- quickloadGenomeId(genome)
  qlg <- new("QuickloadGenome", quickload = quickload, genome = genome_id)
  if (create) {
    if (is.character(genome) && missing(seqinfo))
      stop("No seqinfo for genome '", genome, "'")
    createQuickloadGenome(qlg, seqinfo, title)
  }
  qlg
}

setMethod("show", "QuickloadGenome", function(object) {
  cat(class(object), "track database\ngenome:", genome(object), "\nquickload:",
      uri(quickload(object)), "\n")
  cat(IRanges:::labeledLine("names", names(object)))
})

genomeFile <- function(x) file.path(uri(x), "mod_chromInfo.txt")
annotsFile <- function(x) file.path(uri(x), "annots.xml")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import of tracks/sequence from Quickload
###

setMethod("track", "QuickloadGenome", function(object, name, ...) {
  emd <- elementMetadata(object)
  if (!name %in% emd$title)
    stop("Track '", name, "' does not exist")
  md <- as.list(emd[emd$title == name,])
  rd <- import(file.path(uri(object), md$name), ...)
  metadata(rd)$quickload <- md
  rd
})

## Since a QuickloadGenome can store only a
## single referenece genome, it is best to treat it as a "slot".

setGeneric("referenceSequence",
           function(x, ...) standardGeneric("referenceSequence"))

referenceSequenceFile <- function(x) {
  paste(file.path(uri(x), genome(x)), ".2bit", sep = "")
}

setMethod("referenceSequence", "QuickloadGenome",
          function(x, which = as(seqinfo(x), "GenomicRanges"), ...)
          {
            import(referenceSequenceFile(x), which = which, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export of tracks/sequence to Quickload
###

### FIXME: check for file URI scheme

haveRsamtools <- function() {
  isTRUE(try(packageVersion("Rsamtools") >= "1.7.0", silent = TRUE))
}

indexTrack <- function(track, uri, format) {
  indexed <- NULL
  if (haveRsamtools()) { ## a bunch of heuristics to get an index, if possible
    formats <- eval(formals(Rsamtools::indexTabix)$format)
    parsed_uri <- parseURI(uri)
    multi <- (is(track, "RangedDataList") ||
              is(track, "SimpleGenomicRangesList")) && length(track) > 1
    if (format %in% formats && uriIsLocal(parsed_uri) && !multi) {
      original_path <- parsed_uri$path
      path <- Rsamtools::bgzip(original_path, overwrite = TRUE)
      skip <- if (!is.null(scanTrackLine(uri))) 1L else 0L
      Rsamtools::indexTabix(path, format, skip = skip)
      indexed <- Rsamtools::TabixFile(path)
      unlink(original_path)
    }
  }
  indexed
}

.exportToQuickload <-function(object, name,
                              format = bestFileFormat(value, object),
                              index = TRUE, metadata = character(), ..., value)
{
  if (is(value, "RsamtoolsFile")) {
    if (missing(name))
      name <- basename(path(value))
    track(object, name, metadata = metadata) <- path(value)
    copyResourceToQuickload(object, Rsamtools::index(value))
  } else {
    if (is(object, "Annotated")) {
      value_metadata <- metadata(value)$quickload
      value_metadata[names(metadata)] <- metadata
      metadata <- value_metadata
    }
    filename <- paste(name, format, sep = ".")
    path <- file.path(uri(object), filename)
    seqinfo(value) <- seqinfo(object)
    export(value, path, format = format, ...)
    if (index && !is.null(indexed <- indexTrack(value, path, format)))
      track(object, name, index = FALSE, metadata = metadata) <- indexed
    else track(object, name, metadata = metadata) <- path
  }
  object
}

setReplaceMethod("track",
                 signature(object = "QuickloadGenome", value = "ANY"),
                 .exportToQuickload)

copyResourceToQuickload <- function(object, path) {
  uri <- parseURI(path)
  if (uri$scheme == "")
    path <- paste("file://", path, sep = "")
  filename <- basename(path)
  object_uri <- parseURI(uri(object))
  if (uriIsLocal(object_uri)) {
    dest_file <- file.path(object_uri$path, filename)
    if (dest_file != path)
      download.file(path, dest_file)
  } else stop("Quickload is not local; cannot copy track")
  filename
}

setReplaceMethod("track",
                 signature(object = "QuickloadGenome",
                           value = "character"),
                 function(object, name = basename(object),
                          metadata = character(), value)
                 {
                   file <- URLencode(copyResourceToQuickload(object, value))
                   files <- QuickloadGenome_annotFiles(object)
                   ## make sure we have the three required attributes
                   attrs <- c(name = file, title = name)
                   attrs[names(metadata)] <- metadata
                   filenames <- getNodeSet(files, "//@name")
                   if (file %in% filenames) {
                     removeChildren(files, match(file, filenames))
                   }
                   files <- addChildren(files,
                                        newXMLNode("file", attrs = attrs))
                   saveXML(files, annotsFile(object))
                   object
                 })

setGeneric("referenceSequence<-",
           function(x, ..., value) standardGeneric("referenceSequence<-"))

setReplaceMethod("referenceSequence",
                 signature(x = "QuickloadGenome", value = "ANY"),
                 function(x, name, ..., value)
                 {
                   export.2bit(value, referenceSequenceFile(x))
                   x
                 })

### Not exported yet
setGeneric("synonyms", function(x, ...) standardGeneric("synonyms"))

synonymFile <- function(x) file.path(uri(x), "synonyms.txt")

setMethod("synonyms", "QuickloadGenome", function(x) {
  CharacterList(strsplit(readLines(synonymFile(x)), "\t"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

normArgGenome <- function(genome) {
  if (is.character(genome)) {
    bs_genome <- BSGenomeForID(genome)
    if (!is.null(bs_genome))
      genome <- bs_genome
  }
  genome
}

quickloadGenomeId <- function(genome) {
  if (!is.character(genome)) {
    species <- sub("^(.).*? ", "\\1_", organism(genome))
    date <- sub("\\. ", "_", releaseDate(genome))
    paste(species, date, sep = "_")
  } else genome
}

createQuickloadGenome <- function(x, seqinfo, title) {
  if (genome(x) %in% genome(quickload(x))) {
    message("NOTE: Genome '", genome(x), "' already exists")
    return()
  }
  createResource(uri(x), dir = TRUE)
  createResource(annotsFile(x), content = "<files/>")
  seqinfo(x) <- seqinfo
  addGenomeToContents(x, title)
}
