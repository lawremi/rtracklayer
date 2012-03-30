test_twoBit <- function() {
  test_path <- system.file("tests", package = "rtracklayer")

  test_2bit <- file.path(test_path, "test.2bit")

  correct_name <- "gi|157704452|ref|AC_000143.1| Homo sapiens chromosome 11, alternate assembly (based on HuRef), whole genome shotgun sequence"
  correct_seq <- "TGATGGAAGAATTATTTGAAAGCCATATAGAATGAAATGACTCTATACCCAAATTAAAACTCAAAAACTTACTCAAAATAGTCCAGAGACTACAACTTCA"
  correct_char <- setNames(correct_seq, correct_name)
  correct_2bit <- Biostrings::DNAStringSet(correct_char)

  ## TEST: basic import
  test <- import(test_2bit)
  checkIdentical(test, correct_2bit)

  ## TEST: basic export, import
  test_2bit_out <- file.path(tempdir(), "test_out.2bit")
  export(correct_2bit, test_2bit_out)
  on.exit(unlink(test_2bit_out))
  test <- import(test_2bit_out)
  checkIdentical(test, correct_2bit)

  ## TEST: twoBit extension
  test_twoBit_out <- file.path(tempdir(), "test_out.twoBit")
  export(correct_2bit, test_twoBit_out)
  on.exit(unlink(test_twoBit_out))
  test <- import(test_twoBit_out)
  checkIdentical(test, correct_2bit)

  ## TEST: character export
  export(correct_char, test_2bit_out)
  test <- import(test_2bit_out)
  checkIdentical(test, correct_2bit)
  
  ## TEST: 'which'
  which_range <- IRanges(c(10, 40), c(30, 42))
  correct_which <- seqselect(correct_2bit[[1]], which_range)
  which <- GRanges(names(correct_2bit), which_range)
  test <- import(test_2bit_out, which = which)
  checkIdentical(unlist(test), correct_which)

  ## TEST: empty which
  which_range <- IRanges()
  correct_which <- Biostrings::DNAStringSet()
  which <- GRanges(character(), which_range)
  test <- import(test_2bit_out, which = which)
  checkIdentical(as.character(test), as.character(correct_which))

  ## TEST: which with empty range
  which_range <- IRanges(1, 0)
  which <- GRanges(names(correct_2bit), which_range)
  test <- import(test_2bit_out, which = which)
  checkIdentical(unlist(correct_which), unlist(test))
}
