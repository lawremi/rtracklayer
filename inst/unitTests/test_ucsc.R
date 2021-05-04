test_ucsc <- function(x) {

    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ### TEST ucscTableQuery Class
    ###

    genome <- "hg38"
    table_name <- "gold"
    full_range <- as(Seqinfo(genome = genome), "GRanges")
    custom_range <- GRangesForUCSCGenome(genome, "chr1", IRanges(67003232, 67132477))
    selected_table <- data.frame(bin = 963, chrom = "chr5", chromStart = 49656261,
                                 chromEnd = 49661871, ix = 804, type ="W", frag="ABBA01004242.1",
                                 fragStart = 0, fragEnd = 5610, strand = "+")

    # creating a track and a table for UCSCSession and genome identifier
    elementMetadata  <- list(bin = 0, ix = 1060, type = "F", frag = "AL133320.8",
                             fragStart = 2000, fragEnd = 131245)
    track <- GRanges("chr1", IRanges(67003232, 67132477), "+", elementMetadata)
    genome(track) <- genome
    table <- data.frame(bin = 0, chrom = "chr1", chromStart = 67003232, chromEnd = 67132477,
                        ix = 1060, type = "F", frag = "AL133320.8", fragStart = 2000, 
                        fragEnd = 131245,strand = "+")


    test_trackhub_path <- system.file("tests", "trackhub", package = "rtracklayer")
    trackhub_genome <- "hg19"
    trackhub_custom_range <- GRangesForUCSCGenome("hg19", "chr1", IRanges(237640, 237791))
    trackhub_full_range <- as(Seqinfo(genome = trackhub_genome), "GRanges")
    trackhub_table_name <- "wgEncodeUWDukeDnaseGM12878FdrPeaks"

    # creating a track for Track Hub
    start <- c(237640, 521500 ,565725, 565900, 566760,
               119905, 122525, 173925, 179865, 180185)
    ir <- IRanges(start, width = 151)
    space <- factor(c(rep("chr1", 5), rep("chr10", 5)))
    name <- rep(".", 10)
    score <- seq.int(70L, 700L, length = 10)
    signalValue <- seq(10, 100, length = 10)
    peak <- rep(-1L, 10)
    trackhub_track <- GRanges(space, ir, name = name, score = score,
                              signalValue = signalValue , peak = peak)
    si <- Seqinfo(genome = "hg19")
    seqlengths(trackhub_track) <- seqlengths(si)[levels(space)]

    ## TEST: ucscTableQuery with UCSCSession with NAMES selection
    session <- browserSession()
    query <- ucscTableQuery(session, table = table_name, names = "ABBA01004242.1")
    checkIdentical(range(query), full_range)
    checkIdentical(getTable(query), selected_table)

    ## TEST: ucscTableQuery with UCSCSession without any custom range selection
    query <- ucscTableQuery(session, table = table_name)
    checkIdentical(range(query), full_range)

    ## TEST: ucscTableQuery with UCSCSession with custom range selection
    query <- ucscTableQuery(session, table = table_name, range = custom_range)
    checkIdentical(range(query), custom_range)

    ## TEST: ucscTableQuery with genome idenetifer(character) without any custom range selection
    query <- ucscTableQuery(genome, table = table_name)
    checkIdentical(genome(query), genome)
    checkIdentical(tableName(query), table_name)
    checkIdentical(range(query), full_range)

    ## TEST: ucscTableQuery with genome idenetifer(character) with custom range selection
    query <- ucscTableQuery(genome, table = table_name, range = custom_range)
    checkIdentical(genome(query), genome)
    checkIdentical(tableName(query), table_name)
    checkIdentical(range(query), custom_range)
    checkIdentical(track(query), track)
    checkIdentical(getTable(query), table)

    if (.Platform$OS.type == "windows")
        return()

    ## TEST: ucscTableQuery with Track Hub(character) without any custom range selection
    query <- ucscTableQuery(test_trackhub_path, table = trackhub_table_name,
                            genome = trackhub_genome)
    checkIdentical(genome(query), trackhub_genome)
    checkIdentical(tableName(query), trackhub_table_name)
    checkIdentical(range(query), trackhub_full_range)
    checkIdentical(track(query), trackhub_track)
    checkIdentical(getTable(query), as.data.frame(trackhub_track))

    ## TEST: ucscTableQuery with Track Hub(character) with custom range selection
    query <- ucscTableQuery(test_trackhub_path, table = trackhub_table_name,
                            genome = trackhub_genome, range = trackhub_custom_range)
    checkIdentical(genome(query), trackhub_genome)
    checkIdentical(tableName(query), trackhub_table_name)
    checkIdentical(range(query), trackhub_custom_range)
    checkIdentical(track(query), trackhub_track[1])
    checkIdentical(getTable(query), as.data.frame(trackhub_track)[1,])
}
