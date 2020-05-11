test_trackhub <- function()    {
    test_trackhub_path <- system.file("tests", "trackhub", package = "rtracklayer")
    th <- TrackHub(test_trackhub_path)

    correct_uri <- file.path("file://", test_trackhub_path)
    correct_genome <- "hg19"
    correct_length <- as.integer(1)

    ## TEST: uri
    checkIdentical(uri(th), correct_uri)

    ## TEST: genome
    checkIdentical(genome(th), correct_genome)

    ## TEST: length
    checkIdentical(length(th), correct_length)
}
