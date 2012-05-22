test_bed <- function() {
  test_path <- system.file("tests", package = "rtracklayer")

  ## Import of classic test.bed with RGB colors
  test_bed <- file.path(test_path, "test.bed")
  
  ## import()

  createCorrectRd <- function(seqinfo) {

    ir <- IRanges(c(127471197, 127472364, 127473531, 127474698, 127475865),
                  width = 1167)
    space <- factor(rep(c("chr7", "chr16"), c(3, 2)), seqlevels(seqinfo))
    blocks <- split(IRanges(c(1, 501, 1068, 1, 668, 1, 1, 1),
                            c(300, 700, 1167, 250, 1167, 1167, 1167, 1167)),
                    rep(seq_len(5), c(3, 2, 1, 1, 1)))
    names(blocks) <- NULL
    correct_rd <- RangedData(ir,
                             name = c("Pos1", "Pos2", "Neg1", "Pos3", "Neg2"),
                             score = c(0, 2, 0, 5, 5),
                             strand = strand(c("+", "+", "-", "+", "-")),
                             itemRgb = c("#FF0000", "#FF0000", "#FF0000",
                               "#FF0000", "#0000FF"),
                             thick = ir, blocks, space = space)
    if (!any(is.na(genome(seqinfo))))
      universe(correct_rd) <- unname(genome(seqinfo)[1])
    metadata(ranges(correct_rd))$seqinfo <- seqinfo
    correct_rd
  }

  createCorrectUCSC <- function(rd) {
    new("UCSCData", rd,
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
  
  correct_rd <- createCorrectRd(Seqinfo(c("chr7", "chr16")))
  correct_ucsc <- createCorrectUCSC(correct_rd)
  
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
  close(test_bed_con)
  
  test_bed_con <- file(test_bed, "r")
  test <- import(test_bed_con, format = "bed")
  checkIdentical(test, correct_ucsc)
  close(test_bed_con)

  test_bed_con <- file(test_bed)
  test <- import(BEDFile(test_bed_con))
  checkIdentical(test, correct_ucsc)
  close(test_bed_con)
  
  test <- import(test_bed, trackLine = FALSE)
  checkIdentical(test, correct_rd)

  correct_gr <- as(correct_ucsc, "GRanges")
  test <- import(test_bed, asRangedData = FALSE)
  checkIdentical(correct_gr, sort(test))

  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectUCSC(createCorrectRd(hg19_seqinfo))
  test <- import(test_bed, genome = "hg19")
  checkIdentical(correct_genome, test)

  subcols <- c("name", "strand", "thick")
  correct_subcols <- correct_ucsc[,subcols]
  test <- import(test_bed, colnames = subcols)
  checkIdentical(correct_subcols, test)
  
  which <- RangesList(chr7 = ranges(correct_rd)[[1]][1:2])
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
  correct_empty <- new("UCSCData", RangedData(),
                       trackLine = correct_ucsc@trackLine)
  seqinfo(correct_empty) <- Seqinfo()
  checkIdentical(test, correct_empty)
  close(test_bed_con)
  
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
  correct_list <- RangedDataList(ItemRGBDemo = correct_ucsc,
                                 ItemRGBDemo2 = correct_ucsc2)
  checkIdentical(correct_list, test)

  export(correct_ucsc, test_bed_out, name = "ItemRGBDemo2")
  test <- import(test_bed_out)
  checkIdentical(correct_ucsc2, test)

  test_bed_url <- paste("file:///", test_bed_out, sep = "")
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
  
  ## RangedDataList
  
  export(correct_list, test_bed_out)
  test <- import(test_bed_out)
  checkIdentical(correct_list, test)

  ## GenomicRangesList
  
  test <- import(test_bed_out, asRangedData = FALSE)
  checkIdentical(as(correct_list, "GenomicRangesList"), test)
  
  ## To/From gzip

  test_bed_gz <- paste(test_bed_out, ".gz", sep = "")
  on.exit(unlink(test_bed_gz))
  export(correct_ucsc, test_bed_gz)
  test <- import(test_bed_gz)
  checkIdentical(correct_ucsc, test)

  export(correct_ucsc2, test_bed_gz, append = TRUE)
  test <- import(test_bed_gz)
  checkIdentical(correct_list, test)
  
  test_bed_gz_url <- paste("file:///", test_bed_gz, sep = "")
  export(correct_ucsc, test_bed_gz_url)
  test <- import(test_bed_gz_url)
  checkIdentical(correct_ucsc, test)
  
  ## To/From tabix

  export(correct_ucsc, test_bed_out, index = TRUE)
  on.exit(unlink(paste(test_bed_gz, ".tbi", sep = "")))
  test <- import(test_bed_gz, which = which)
  checkIdentical(correct_which, test)

  ## check TabixFile
  
  test_bed_tabix <- Rsamtools::TabixFile(test_bed_gz)
  test <- import(test_bed_tabix)
  checkIdentical(correct_ucsc, test)

  ## look mom, no track line
  export(correct_ucsc, test_bed_out, index = TRUE, trackLine = FALSE)
  test <- import(test_bed_gz, which = which)
  checkIdentical(subsetByOverlaps(correct_rd, which), test)
  test <- import(test_bed_tabix, format = "foo")
  
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
    
  ## Set seqinfo on correct_rd, coerce to UCSCData, export, import and check

  correct_genome_rd <- as(correct_genome, "RangedData")
  
  ## Set offset in correct_ucsc track line, export, then:
   
  correct_offset <- correct_rd
  ranges(correct_offset) <- shift(ranges(correct_offset), -1)
  correct_ucsc@trackLine@offset <- 1L
  export(correct_ucsc, test_bed_out)
  test <- import(test_bed_out, trackLine = FALSE)
  checkIdentical(test, correct_offset)

  test <- import(test_bed_out)
  checkIdentical(test, correct_ucsc)
  correct_ucsc@trackLine@offset <- 0L
  
  ## Drop all extra columns, see if it still works
  correct_stripped <- correct_ucsc[,character()]
  export(correct_stripped, test_bed_out)
  test <- import(test_bed_out)
  checkIdentical(test, correct_stripped)
  
  ## - and even when asking for a column like blocks:
  correct_blocks <- correct_ucsc[,"blocks"]
  test <- import(test_bed, colnames = "blocks")
  checkIdentical(test, correct_blocks)
  
  ## Drop the columns except for blocks, see if it still works
  export(correct_blocks, test_bed_out)
  test <- import(test_bed_out)
  correct_fill_to_blocks <- correct_ucsc
  correct_fill_to_blocks <- within(correct_fill_to_blocks, {
    name <- NA_character_
    score <- 0
    strand <- strand("*")
    itemRgb <- NA_character_
    thick <- unlist(ranges, use.names = FALSE)
  })
  checkIdentical(test, correct_fill_to_blocks)
}
