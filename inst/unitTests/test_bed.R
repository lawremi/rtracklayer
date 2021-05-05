test_bed <- function() {
  test_path <- system.file("tests", package = "rtracklayer")

  ## Import of classic test.bed with RGB colors
  test_bed <- file.path(test_path, "test.bed")
  
  ## import()

  createCorrectGR <- function(seqinfo) {
    ir <- IRanges(c(127471197, 127472364, 127473531, 127474698, 127475865),
                  width = 1167)
    space <- factor(rep(c("chr7", "chr9"), c(3, 2)), seqlevels(seqinfo))
    blocks <- split(IRanges(c(1, 501, 1068, 1, 668, 1, 1, 1),
                            c(300, 700, 1167, 250, 1167, 1167, 1167, 1167)),
                    rep(seq_len(5), c(3, 2, 1, 1, 1)))
    names(blocks) <- NULL
    correct_gr <- GRanges(space, ir,
                          strand = strand(c("+", "+", "-", "+", "-")),
                          name = c("Pos1", "Pos2", "Neg1", "Pos3", "Neg2"),
                          score = c(0, 2, 0, 5, 5),
                          itemRgb = c("#FF0000", "#FF0000", "#FF0000",
                            "#FF0000", "#0000FF"),
                          thick = ir, blocks)
    seqinfo(correct_gr) <- seqinfo
    correct_gr
  }

  createCorrectUCSC <- function(gr) {
    new("UCSCData", gr,
        trackLine = new("BasicTrackLine", itemRgb = TRUE,
          name = "ItemRGBDemo",
          description = "Item RGB demonstration",
          visibility = "2", color = c(0L, 60L, 120L),
          priority = 1, group = "user", offset = 0L,
          url = "http://genome.ucsc.edu/foo.html?query=$$",
          htmlUrl = paste("http://genome.ucsc.edu/goldenPath/",
            "help/ct_description.txt", sep = ""),
          colorByStrand =
          matrix(c(0L, 0L, 255L, 255L, 0L, 0L), 3),
          useScore = FALSE))
  }
  
  correct_gr <- createCorrectGR(Seqinfo(c("chr7", "chr9")))
  correct_ucsc <- createCorrectUCSC(correct_gr)
  
  test <- import(test_bed)
  checkIdentical(test, correct_ucsc)
  
  test_bed_file <- BEDFile(test_bed)
  test <- import(test_bed_file)
  checkIdentical(test, correct_ucsc)
  checkIdentical(import(test_bed_file, format = "bed"), correct_ucsc)
  checkException(import(test_bed_file, format = "gff"))

  test_bed_con <- file(test_bed)
  test <- import(test_bed_con, format = "bed")
  checkIdentical(test, correct_ucsc)
  
  test_bed_con <- file(test_bed, "r")
  test <- import(test_bed_con, format = "bed")
  checkIdentical(test, correct_ucsc)
  close(test_bed_con)

  test_bed_con <- file(test_bed)
  test <- import(BEDFile(test_bed_con))
  checkIdentical(test, correct_ucsc)
  
  test <- import(test_bed, trackLine = FALSE)
  checkIdentical(test, correct_gr)

  test <- import(test_bed)
  checkIdentical(correct_ucsc, test)

  if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
    stop("'BSgenome.Hsapiens.UCSC.hg19' must be installed to run tests")
  }
  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectUCSC(createCorrectGR(hg19_seqinfo))
  test <- import(test_bed, genome = "hg19")
  checkIdentical(correct_genome, test)

  subcols <- c("name", "thick")
  correct_subcols <- correct_ucsc
  mcols(correct_subcols) <- mcols(correct_subcols)[ , subcols]
  test <- import(test_bed, colnames = c(subcols, "strand"))
  checkIdentical(correct_subcols, test)
  strand(correct_subcols) <- "*"
  test <- import(test_bed, colnames = subcols)
  checkIdentical(correct_subcols, test)
  
  which <- correct_gr[1:2]
  correct_which <- subsetByOverlaps(correct_ucsc, which)
  test <- import(test_bed, which = which)
  checkIdentical(correct_which, test)

  test <- import(test_bed, format = "bed")
  checkIdentical(correct_ucsc, test)
  
  ## import.bed()

  test <- import.bed(test_bed)
  checkIdentical(correct_ucsc, test)

  test_bed_con <- pipe(paste("head -n2", test_bed))
  test <- import.bed(test_bed_con)
  correct_empty <- new("UCSCData", GRanges(),
                       trackLine = correct_ucsc@trackLine)
  seqinfo(correct_empty) <- Seqinfo()
  checkIdentical(test, correct_empty)
  
  ## export()

  ## the 'gsub' is to handle Windows paths (for later coercion to URL)
  test_bed_out <- gsub("\\\\", "/", file.path(tempdir(), "test.bed"))
  on.exit(unlink(test_bed_out))
  export(correct_ucsc, test_bed_out)
  test <- import(test_bed_out)
  checkIdentical(correct_ucsc, test)

  export(correct_ucsc, test_bed_out, format = "bed")
  test <- import(test_bed_out)
  checkIdentical(correct_ucsc, test)

  test_bed_out_file <- BEDFile(test_bed_out)
  export(correct_ucsc, test_bed_out_file)
  test <- import(test_bed_out)
  checkIdentical(correct_ucsc, test)
  checkException(export(correct_ucsc, test_bed_out_file, format = "gff"))

  correct_ucsc2 <- initialize(correct_ucsc,
                              trackLine = initialize(correct_ucsc@trackLine,
                                                     name = "ItemRGBDemo2"))
  export(correct_ucsc2, test_bed_out_file, append = TRUE)
  test <- import(test_bed_out_file)
  correct_list <- GRangesList(ItemRGBDemo = correct_ucsc,
                              ItemRGBDemo2 = correct_ucsc2,
                              compress=FALSE)
  checkIdentical(correct_list, test)

  export(correct_ucsc, test_bed_out, name = "ItemRGBDemo2")
  test <- import(test_bed_out)
  checkIdentical(correct_ucsc2, test)

  test_bed_url <- paste("file://", test_bed_out, sep = "")
  export(correct_ucsc, test_bed_url)
  test <- import(test_bed_url)
  checkIdentical(correct_ucsc, test)

if (FALSE) { # enable to test an HTTP URL using the R help server
  http_pipe <- pipe("Rscript -e 'tools::startDynamicHelp(); writeLines(as.character(tools:::httpdPort)); Sys.sleep(10);' &")
  port <- readLines(http_pipe, n = 1)
  test_bed_http <- paste("http://127.0.0.1:", port,
                         "/library/rtracklayer/doc/example.bed", sep = "")
  test <- import(test_bed_http)
  checkIdentical(correct_ucsc, test)
  close(http_pipe)
}
  
  ## GenomicRangesList
  
  export(correct_list, test_bed_out)
  test <- import(test_bed_out)
  checkIdentical(correct_list, test)

  ## To/From gzip

  test_bed_gz <- paste(test_bed_out, ".gz", sep = "")
  on.exit(unlink(test_bed_gz))
  export(correct_ucsc, test_bed_gz)
  test <- import(test_bed_gz)
  checkIdentical(correct_ucsc, test)

  export(correct_ucsc2, test_bed_gz, append = TRUE)
  test <- import(test_bed_gz)
  checkIdentical(correct_list, test)
  
  test_bed_gz_url <- paste("file://", test_bed_gz, sep = "")
  export(correct_ucsc, test_bed_gz_url)
  test <- import(test_bed_gz_url)
  checkIdentical(correct_ucsc, test)
  
  ## To/From tabix

  test_bed_bgz <- paste(test_bed_out, ".bgz", sep = "")
  export(correct_ucsc, test_bed_out, index = TRUE)
  on.exit(unlink(paste(test_bed_bgz, ".tbi", sep = "")))
  test <- import(test_bed_bgz, which = which)
  checkIdentical(correct_which, test)

  ## check TabixFile
  
  test_bed_tabix <- Rsamtools::TabixFile(test_bed_bgz)
  test <- import(test_bed_tabix)
  checkIdentical(correct_ucsc, test)

  ## look mom, no track line
  export(correct_ucsc, test_bed_out, index = TRUE, trackLine = FALSE)
  test <- import(test_bed_bgz, which = which)
  checkIdentical(subsetByOverlaps(correct_gr, which), test)
  #test <- import(test_bed_tabix, format = "foo")
  
  ## To/From text
  
  bed_text <- export(correct_ucsc, format = "bed")
  test <- import(format = "bed", text = bed_text)
  checkIdentical(correct_ucsc, test)

  ## TODO: empty text
  
  ## Using connection to add comment header

  test_bed_con <- file(test_bed_out)
  open(test_bed_con, "w")
  comment <- "# test comment"
  writeLines(comment, test_bed_con)
  export(correct_ucsc, test_bed_con)
  close(test_bed_con)
  checkIdentical(comment, readLines(test_bed_out, n = 1))
  test <- import(test_bed_out)
  checkIdentical(correct_ucsc, test)
    
  ## Set offset in correct_ucsc track line, export, then:
   
  correct_offset <- shift(correct_gr, -1)
  correct_ucsc@trackLine@offset <- 1L
  export(correct_ucsc, test_bed_out)
  test <- import(test_bed_out, trackLine = FALSE)
  checkIdentical(test, correct_offset)

  test <- import(test_bed_out)
  checkIdentical(test, correct_ucsc)
  correct_ucsc@trackLine@offset <- 0L
  
  ## Drop all extra columns, see if it still works
  correct_stripped <- correct_ucsc
  mcols(correct_stripped) <- NULL
  export(correct_stripped, test_bed_out)
  test <- import(test_bed_out)
  mcols(correct_stripped)$name <- NA_character_
  mcols(correct_stripped)$score <- 0
  checkIdentical(test, correct_stripped)
  
  ## - and even when asking for a column like blocks:
  correct_blocks <- correct_ucsc
  mcols(correct_blocks) <- mcols(correct_blocks)[ , "blocks", drop=FALSE]
  test <- import(test_bed, colnames = c("blocks", "strand"))
  checkIdentical(test, correct_blocks)
  strand(correct_blocks) <- "*"
  test <- import(test_bed, colnames = "blocks")
  checkIdentical(test, correct_blocks)
  
  ## Drop the columns except for blocks, see if it still works
  export(correct_blocks, test_bed_out)
  test <- import(test_bed_out)
  correct_fill_to_blocks <- correct_ucsc
  strand(correct_fill_to_blocks) <- "*"
  mcols(correct_fill_to_blocks)$name <- NA_character_
  mcols(correct_fill_to_blocks)$score <- 0
  mcols(correct_fill_to_blocks)$itemRgb <- NA_character_
  mcols(correct_fill_to_blocks)$thick <- ranges(correct_fill_to_blocks)
  checkIdentical(test, correct_fill_to_blocks)
}

test_extendedBed <- function()
{

        # the narrowPeak format represents a variety of "BED6+N" formats
        # used by the ENCODE project.  see
        # http://genome.ucsc.edu/FAQ/FAQformat.html
        # "This format is used to provide called peaks of signal
        #  enrichment based on pooled, normalized (interpreted) data.
        #  It is a BED6+4 format."
  
    file <- system.file("extdata", "demo.narrowPeak.gz",  package="rtracklayer")
    extraCols <- c(signalValue="numeric", pValue="numeric", qValue="numeric",
                   peak="integer")
    gr <- import(file, forma="bed", extraCols=extraCols, genome="hg19")
    checkEquals(length(gr), 6)
    checkEquals(colnames(mcols(gr)),
                c("name","score","signalValue","pValue","qValue","peak"))

        # make sure that all seqnames in the gr object
        # are also in the seqinfo(gr) object
    checkTrue(all(seqnames(gr) %in% names(seqinfo(gr))))


    
} 

test_bedpe <- function() {
    path <- system.file("tests", "test.bedpe", package="rtracklayer")
    nms <- c("TUPAC_0001:3:1:0:1452#0", "TUPAC_0001:3:1:0:1472#0",
             "TUPAC_0001:3:1:1:1833#0")
    gr1 <- GRanges(c("chr7", "chr11", "chr20"),
                   IRanges(c(118965073, 46765607, 54704675),
                           c(118965122, 46765656, 54704724)),
                   strand="+")
    gr2 <- GRanges(c("chr7", "chr10", "chr20"),
                   IRanges(c(118970080, 46769935, 54708988),
                           c(118970129, 46769984, 54709037)),
                   strand="-")
    seqlevels(gr1) <- union(seqlevels(gr1), seqlevels(gr2))
    seqlevels(gr2) <- seqlevels(gr1)
    pairs <- Pairs(gr1, gr2, name=nms, score=37)
    bedpe <- import(path)

    checkIdentical(pairs, bedpe)

    # test export
    test_bedpe_out <- file.path(tempdir(), "test.bedpe")
    on.exit(unlink(test_bedpe_out))
    export(bedpe, test_bedpe_out)
    test <- import(test_bedpe_out)
    checkIdentical(bedpe, test)

    bedpe2 <- bedpe
    mcols(bedpe2)$qvalue <- c(0.02, 0.03, 0.05)
    mcols(bedpe2)$annotation <- c("promoter", "enhancer", "intron")

    # test extended bedpe
    test_bedpe2_out <- file.path(tempdir(), "test2.bedpe")
    on.exit(unlink(test_bedpe2_out))
    export(bedpe2, test_bedpe2_out)
    test <- import(test_bedpe2_out,
                   extraCols=c(qvalue="numeric", annotation="character"))
    checkIdentical(bedpe2, test)
}
