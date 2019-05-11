test_gff <- function() {
  test_path <- system.file("tests", package = "rtracklayer")

  space <- c(rep("chr10", 15), rep("chr12", 16))
  start <- c(rep(92828, 4), 92997, rep(94555, 4), rep(94744, 3), rep(95122, 2),
             95348, rep(87984, 4),
             rep(c(88257, 88570, 88860, 89675, 90587, 90796), each = 2))
  end <- c(95504, 95178, 95504, rep(94054, 2), rep(94665, 2), 94615, 94665,
           rep(94852, 3), rep(95178, 2), 95504,
           rep(c(91263, 88017, 88392, 88771, 89018, 89827, 90655, 91263),
               each = 2))
  type <- c("gene", "mRNA", "mRNA", "exon", "CDS", "exon", "exon",
            "CDS", "CDS", "exon", "exon", "CDS", "exon", "CDS", "exon",
            "gene", "mRNA", "exon", "CDS", "exon", "CDS", "exon", "CDS",
            "exon", "CDS", "exon", "CDS", "exon", "CDS", "exon", "CDS")
  type <- factor(type, unique(type))
  source <- factor("rtracklayer")
  phase <- NA_integer_
  score <- c(5, rep(NA, length(type) - 1L))
  strand <- strand(c(rep("-", 14), "*", rep("+", 15), "*"))
  Alias <- CharacterList(c(list(c("FLJ40100", "TUBB8")),
                           rep(list(character()), 14), "LOC100288778",
                           rep(list(character()), 15)))
  ID <- c("GeneID:347688", "873", "872", rep(NA, 12), "GeneID:100288778",
          "4644", rep(NA, 14))
  Name <- c(rep("TUBB8", 3), rep(NA, 12), rep("LOC100288778", 2),
            rep(NA, 14))
  Parent <- CharacterList(c(list(character()), rep("GeneID:347688", 2),
                            rep(list(c("872", "873")), 2),
                            rep(c("872", "873"), 3), rep("873", 3), "872",
                            list(character()), "GeneID:100288778",
                            rep("4644", 14)))
  geneName <- c("tubulin, beta 8", rep(NA, 14),
                "WAS protein family homolog 1; pseudogene", rep(NA, 15))
  genome <- c("hg19", rep(NA, length(geneName) - 1))
  correct_gff3 <- GRanges(space, IRanges(start, end), strand,
                          source, type, score, phase,
                          ID, Name, geneName, Alias, genome, Parent)
  seqinfo(correct_gff3) <- Seqinfo(c("chr10", "chr12"))

  correct_gff1 <- correct_gff3[,c("source", "type", "score", "phase")]
  correct_gff1$group <- as.factor(seqnames(correct_gff3))

  correct_gff2 <- correct_gff3
  toCSV <- function(x) {
    csv <- sapply(x, paste, collapse = ",")
    csv[nchar(csv) == 0] <- NA
    csv
  }
  correct_gff2$Alias <- toCSV(correct_gff2$Alias)
  correct_gff2$Parent <- toCSV(correct_gff2$Parent)

  ## TEST: basic GFF3 import
  test_gff3 <- file.path(test_path, "genes.gff3")
  test <- import(test_gff3)
  checkIdentical(correct_gff3, test)

  ## TEST: import.gff*
  test <- import.gff(test_gff3)
  checkIdentical(correct_gff3, test)
  test <- import.gff3(test_gff3)
  checkIdentical(correct_gff3, test)
  #suppressWarnings(test <- import.gff2(test_gff3))
  #checkIdentical(correct_gff3, test)
  oldOpts <- options(warn = 2)
  checkException(import.gff2(test_gff3))
  options(oldOpts)
  
  ## TEST: GFF(3)File
  test_gff_file <- GFF3File(test_gff3)
  test <- import(test_gff_file)
  checkIdentical(correct_gff3, test)
  test_gff_file <- GFFFile(test_gff3)
  test <- import(test_gff_file)
  checkIdentical(correct_gff3, test)
  test_gff_file <- GFFFile(test_gff3, version = "3")
  test <- import(test_gff_file)
  checkIdentical(correct_gff3, test)
  test_gff_file <- GFF2File(test_gff3)
  #suppressWarnings(test <- import(test_gff_file))
  #checkIdentical(correct_gff3, test)
  oldOpts <- options(warn = 2)
  checkException(test <- import(test_gff_file))
  options(oldOpts)

  ## TEST: 'gff' extension
  test_gff_out <- file.path(tempdir(), "genes.gff")
  on.exit(unlink(test_gff_out))
  export(correct_gff3, test_gff_out)
  test <- import(test_gff_out)
  checkIdentical(test, correct_gff1)
  export(correct_gff3, test_gff_out, version = "1")
  test <- import(test_gff_out)
  checkIdentical(test, correct_gff1)
  export(correct_gff3, test_gff_out, version = "2")
  test <- import(test_gff_out)
  checkIdentical(test, correct_gff2)
  export(correct_gff3, test_gff_out, version = "3")
  test <- import(test_gff_out)
  checkIdentical(test, correct_gff3)
  test <- import(GFF3File(test_gff_out))
  checkIdentical(test, correct_gff3)
  test <- import(GFFFile(test_gff_out))
  checkIdentical(test, correct_gff3)
  test <- import(test_gff_out, version = "3")
  checkIdentical(test, correct_gff3)
  #suppressWarnings(test <- import(test_gff_out, version = "2"))
  #checkIdentical(test, correct_gff3)
  oldOpts <- options(warn = 2)
  checkException(test <- import(test_gff_out, version = "2"))
  options(oldOpts)

  ## TEST: 'gff2' extension
  test_gff2_out <- file.path(tempdir(), "genes.gff2")
  export(correct_gff3, test_gff2_out)
  test <- import(test_gff2_out)
  checkIdentical(test, correct_gff2)

  ## TEST: 'gff1' extension
  test_gff1_out <- file.path(tempdir(), "genes.gff1")
  export(correct_gff3, test_gff1_out)
  test <- import(test_gff1_out)
  checkIdentical(test, correct_gff1)

  ## TEST: 'format' argument
  test_gff_file <- GFF3File(test_gff3)
  test <- import(test_gff_file, format = "gff")
  checkIdentical(test, correct_gff3)
  test <- import(test_gff_file, format = "gff3")
  checkIdentical(test, correct_gff3)
  checkException(import(test_gff_file, format = "gff2"))
  checkException(import(test_gff_file, format = "bed"))
  
  ## TEST: 'genome'  
  si_hg19 <- SeqinfoForBSGenome("hg19")
  correct_hg19 <- correct_gff3
  seqlevels(correct_hg19) <- seqlevels(si_hg19)
  seqinfo(correct_hg19) <- si_hg19
  test <- import(test_gff3, genome = "hg19")
  checkIdentical(test, correct_hg19)
  
  test_gff3_out <- file.path(tempdir(), "genes.gff3")
  on.exit(unlink(test_gff3_out))
  correct_genome_hg19 <- correct_gff3
  genome(correct_genome_hg19) <- "hg19"
  correct_genome_hg19 <- as(correct_genome_hg19, "GRanges")
  export(correct_genome_hg19, test_gff3_out)
  test <- import(test_gff3_out)
  checkIdentical(test, correct_hg19)

  ## TEST: colnames empty, colnames := "geneName"
  test <- import(test_gff3, colnames = character())
  target <- correct_gff3[,character()]
  checkIdentical(target, test)
  test <- import(test_gff3, colnames = "geneName")
  target <- correct_gff3[,"geneName"]
  checkIdentical(target, test)

  ## TEST: import from connection
  test_gff_con <- file(test_gff_out)
  test <- import(test_gff_con, format = "gff")
  checkIdentical(correct_gff3, test)

  ## TEST: export to connection, with preceding comment
  test_gff_con <- file(test_gff_out)
  open(test_gff_con, "w")
  comment <- "# test comment"
  writeLines(comment, test_gff_con)
  export(correct_gff3, test_gff_con, version = "3")
  close(test_gff_con)
  checkIdentical(comment, readLines(test_gff_out, n = 1))
  test <- import(test_gff_out)
  checkIdentical(correct_gff3, test)
  
  ## TEST: 'append'
  export(correct_gff3[seqnames(correct_gff3) == "chr10", ], test_gff3_out)
  export(correct_gff3[seqnames(correct_gff3) == "chr12", ], test_gff3_out,
         append = TRUE)
  test <- import(test_gff3_out)
  checkIdentical(correct_gff3, test)

  ## TEST: 'source'
  target <- correct_gff3
  mcols(target)$source <- factor("test")
  export(correct_gff3, test_gff3_out, source = "test")
  test <- import(test_gff3_out)
  checkIdentical(target, test)

  ## TEST: 'which'
  which <- GRanges("chr10:90000-93000")
  which_target <- subsetByOverlaps(correct_gff3, which)
  test <- import(test_gff3, which = which)
  checkIdentical(which_target, test)
  
  ## TEST: 'index'
  export(correct_gff3, test_gff3_out, index = TRUE)
  test_gff_bgz <- paste(test_gff3_out, ".bgz", sep = "")
  on.exit(unlink(test_gff_bgz))
  on.exit(unlink(paste(test_gff_bgz, ".tbi", sep = "")))
  test <- import(test_gff_bgz, which = which)
  checkIdentical(which_target, test)

  ## TEST: SimpleGRangesList
  ucsc_data1 <- new("UCSCData", keepSeqlevels(correct_gff3, "chr10",
                                              pruning.mode="coarse"),
                    trackLine = new("BasicTrackLine", name = "chr10"))
  ucsc_data2 <- new("UCSCData", keepSeqlevels(correct_gff3, "chr12",
                                              pruning.mode="coarse"),
                    trackLine = new("BasicTrackLine", name = "chr12"))
  correct_grl <- GRangesList(ucsc_data1, ucsc_data2, compress=FALSE)
  mcols(correct_grl[[2]])$genome <- NULL
  names(correct_grl) <- seqlevels(correct_gff3)
  export(correct_grl, test_gff3_out)
  test <- import.ucsc(test_gff3_out)
  checkIdentical(correct_grl, test)
}
