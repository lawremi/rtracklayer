test_trackhub <- function() {
    test_trackhub_path <- system.file("tests", "trackhub", package = "rtracklayer")

    correct_trackhub_uri <- test_trackhub_path
    correct_trackhub_genome <- "hg19"
    correct_trackhub_length <- 1L
    correct_hub <- "test_hub"
    correct_shortLabel <- "test_hub"
    correct_longLabel <- "test_hub"
    correct_genomesFile <- "genomes.txt"
    correct_email <- "user@domain.com"
    correct_descriptionUrl <- "http://www.somedomain.com/articles/h19"
    correct_trackDb <- "hg19/trackDb.txt"

    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ### TEST TrackHub Class
    ###

    th <- TrackHub(test_trackhub_path)

    ## TEST: uri
    checkIdentical(uri(th), correct_trackhub_uri)

    ## TEST: genome
    checkIdentical(genome(th), correct_trackhub_genome)

    ## TEST: length
    checkIdentical(length(th), correct_trackhub_length)

    # TEST: hub
    checkIdentical(hub(th), correct_hub)

    # TEST: shortLabel
    checkIdentical(shortLabel(th), correct_shortLabel)

    # TEST: longLabel
    checkIdentical(longLabel(th), correct_longLabel)

    # TEST: genomesFile
    checkIdentical(genomesFile(th), correct_genomesFile)

    # TEST: email
    checkIdentical(email(th), correct_email)

    # TEST: descriptionUrl
    checkIdentical(descriptionUrl(th), correct_descriptionUrl)

    # TEST: hub<-
    new_hub <- "new_hub"
    hub(th) <- new_hub
    checkIdentical(hub(th), new_hub)
    hub(th) <- correct_hub

    # TEST: shortLabel<-
    new_shortLabel <- "new_hub"
    shortLabel(th) <- new_shortLabel
    checkIdentical(shortLabel(th), new_shortLabel)
    shortLabel(th) <- correct_shortLabel

    # TEST: longLabel<-
    new_longLabel <- "new_hub"
    longLabel(th) <- new_longLabel
    checkIdentical(longLabel(th), new_longLabel)
    longLabel(th) <- correct_longLabel

    # TEST: genomesFile<-
    new_genomesFile <- "newfile.txt"
    genomesFile(th) <- new_genomesFile
    checkIdentical(genomesFile(th), new_genomesFile)
    genomesFile(th) <- correct_genomesFile

    # TEST: email<-
    new_email <- "new@domail.com"
    email(th) <- new_email
    checkIdentical(email(th), new_email)
    email(th) <- correct_email

    # TEST: descriptionUrl<-
    new_descriptionUrl <- "http://newdomail.com/articles/hg19"
    descriptionUrl(th) <- new_descriptionUrl
    checkIdentical(descriptionUrl(th), new_descriptionUrl)
    descriptionUrl(th) <- correct_descriptionUrl

    # TEST: genomeField
    checkIdentical(genomeField(th, "hg19", "trackDb"), correct_trackDb)

    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ### TEST TrackHubGenome Class
    ###

    correct_trackhubgenome_uri <- paste0(correct_trackhub_uri, "/hg19")
    correct_trackhubgenome_genome_name <- "hg19"
    correct_trackhubgenome_length <- 1L
    correct_trackhubgenome_organism <- "BigFoot"
    correct_trackhubgenome_names <- "wgEncodeUWDukeDnaseGM12878FdrPeaks"
    correct_bigDataUrl <- "hg19/wgEncodeUWDukeDnaseGM12878.fdr01peaks.hg19.bb"

    thg <- TrackHubGenome(th, "hg19")

    # TEST: uri
    checkIdentical(uri(thg), correct_trackhubgenome_uri)

    # TEST: genome
    checkIdentical(genome(thg), correct_trackhubgenome_genome_name)

    # TEST: length
    checkIdentical(length(thg), correct_trackhubgenome_length)

    # TEST: organism
    checkIdentical(organism(thg), correct_trackhubgenome_organism)

    # TEST: names
    checkIdentical(names(thg), correct_trackhubgenome_names)

    # TEST: trackNames
    checkIdentical(trackNames(thg), correct_trackhubgenome_names)

    # TEST: trackField
    checkIdentical(trackField(thg, "wgEncodeUWDukeDnaseGM12878FdrPeaks", "bigDataUrl"), correct_bigDataUrl)

    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ### TEST TrackContainer Class
    ###

    correct_slot_type <- "Track"
    correct_track <- Track(track = "tcell", bigDataUrl = "tcell/data.bigWig")

    # TEST: slot type
    tc <- TrackContainer()
    slot_type <- slot(tc, "elementType")
    checkIdentical(slot_type, correct_slot_type)

    # TEST: wrong type slot error reporting
    checkException(tc[[1]] <- 1)

    # TEST: names()
    tc[[1]] <- correct_track
    checkIdentical(names(tc), correct_track@track)
}
