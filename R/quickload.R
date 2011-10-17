### =========================================================================
### Quickload support
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Classes representing Quickload resources
###

setClass("Quickload", representation(uri = "character"), contains = "TrackDb")

uri <- function(x, ...) x@uri

.Quickload_contents <- function(x) {
  contents_file <- file.path(uri(x), "contents.txt")
  if (file.exists(contents_file))
    read.table(contents_file, sep = "\t")
  else NULL
}

setMethod("genome", "Quickload", function(x) {
  contents <- .Quickload_contents(x)
  structure(contents[[1]], names = contents[[2]])
})

Quickload <- function(uri) {
  new("Quickload", uri = normURI(uri))
}

setAs("character", "Quickload", function(from) Quickload(from))

setClass("QuickloadGenome",
         representation(quickload = "Quickload",
                        genome = "character"))

setMethod("uri", "QuickloadGenome",
          function(x) file.path(uri(quickload(x)), genome(x)))

quickload <- function(x, ...) x@quickload

setMethod("genome", "QuickloadGenome", function(x) x@genome)

setMethod("seqinfo", "QuickloadGenome", function(x) {
  chromInfo <- read.table(file.path(uri(x), "mod_chromInfo.txt"), sep = "\t")
  Seqinfo(rownames(chromInfo), chromInfo[[1]], genome = genome(x))
})

setMethod("releaseDate", "QuickloadGenome", function(x) {
  sub(".*?_(.*?)_([^_]*)$", "\\1 \\2", genome(x))
})

setMethod("organism", "QuickloadGenome", function(x) {
  gsub("_", " ", sub("(.*?)_.*?_[^_]*$", "\\1", genome(x)))
})

.QuickloadGenome_annotFiles <- function(x) {
  annots_file <- file.path(uri(x), "annots.xml")
  if (file.exists(annots_file))
    files <- xmlChildren(xmlInternalTreeParse(annots_file))$files
  else files <- xmlNode("files")
}

setMethod("length", "QuickloadGenome", function(x) {
  length(names(x))
})

setMethod("names", "QuickloadGenome", function(x) {
  files <- .QuickloadGenome_annotFiles(x)
  file_names <- getNodeSet(files, "//@name")
  file_titles <- getNodeSet(files, "//@title")
  structure(filenames, names = file_titles)
})

setMethod("elementMetadata", "QuickloadGenome", function(x) {
  files <- .QuickloadGenome_annotFiles(x)
  Reduce(function(x, y) merge(as.data.frame(x), as.data.frame(y), all = TRUE),
         lapply(xmlChildren(files), xmlAttrs))
})

setMethod("track", "QuickloadGenome", function(object, name, ...) {
  tmd <- elementMetadata(object)
  if (!name %in% tmd$name)
    stop("Track '", name, "' does not exist")
  gr <- import(file.path(uri(object), name), ...)
  metadata(gr)$quickload <- as.list(tmd[tmd$name == name,])
  gr
})

QuickloadGenome <- function(quickload, genome) {
  quickload <- as(quickload, "Quickload")
  genome_id <- quickloadGenomeId(genome)
  new("QuickloadGenome", quickload = quickload, genome = genome)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export of genome sequence and range datasets to (local) Quickload
###

setMethod("export", c("BSgenome", "Quickload"),
          function(object, con, igb_jar) {
            genome_id <- quickloadGenomeId(object)
            genome_dir <- quickloadGenomeDir(con, genome_id)
            write.BSgenome(object, genome_dir, "bnib", igb_jar)
            df <- as.data.frame(seqinfo(object))[,1,drop=FALSE]
            write.table(df, file.path(genome_dir, "mod_chromInfo.txt"),
                        quote = FALSE, header = FALSE, sep = "\t")
            genome_dir
          })

setMethod("export", c("GenomicRangesORGRangesList",  "Quickload"),
          function(object, con, file, ...) {
            genome <- singleGenome(genome(object))
            if (length(genome) != 1)
              stop("'object' must be on a single genome")
            genome_id <- quickloadGenomeId(genome)
            genome_dir <- quickloadGenomeDir(con, genome_id)
            data_file <- file.path(genome_dir, file)
            export(object, data_file)
            export(data_file, con, genome_id)
          })

setMethod("export", c("character", "Quickload"),
          function(object, con, genome_id) {
            genome_dir <- quickloadGenomeDir(con, genome_id)
            file <- basename(object)
            dest_file <- file.path(genome_dir, file)
            if (dest_file != object)
              download.file(object, dest_file)
            .QuickloadGenome_annotFiles(genome_dir)
            ## make sure we have the three required attributes
            attrs <- list(name = file, title = file, description = file)
            args <- list(...)
            attrs[names(args)] <- args  
            files <- addChildren(files, xmlNode("file", attrs = attrs))
            saveXML(files, annots_file)
            dest_file
          })

quickloadGenomeDir <- function(quickload, genome_id)
{
  genome_dir <- file.path(uri(quickload), genome_id)
  dir.create(genome_dir, recursive = TRUE)
  
  contents <- .Quickload_contents(quickload)
  if (!genome_id %in% contents[[1]]) {
    desc <- paste(organism(genome), providerVersion(genome))
    contents <- rbind(contents, data.frame(genome_id, desc))
    write.table(contents, contents_file, quote = FALSE, header = FALSE,
                sep = "\t")
  }
  
  genome_dir  
}

quickloadGenomeId <- function(genome) {
  if (is.character(genome))
    genome <- BSGenomeForId(genome)
  if (is(genome, "BSGenome")) {
    species <- sub("BSgenome\\.([^.]*)\\..*", "\\1", bsgenomeName(genome))
    date <- sub("\\. ", "_", releaseDate(genome))
    paste(species, date, sep = "_")
  } else genome
}
