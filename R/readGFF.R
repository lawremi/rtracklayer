### =========================================================================
### readGFF()
### -------------------------------------------------------------------------


.make_filexp_from_filepath <- function(filepath)
{
    if (isSingleString(filepath))
        return(XVector:::open_input_files(filepath)[[1L]])
    if (!is(filepath, "connection"))
        stop(wmsg("'filepath' must be a single string or a connection"))
    if (!isSeekable(filepath))
        stop(wmsg("connection is not seekable"))
    filepath
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGFFPragmas()
###

readGFFPragmas <- function(filepath)
{
    filexp <- .make_filexp_from_filepath(filepath)
    if (is(filexp, "connection")) {
        if (!isOpen(filexp)) {
            open(filexp)
            on.exit(close(filexp))
        }
        if (seek(filexp) != 0) {
            warning(wmsg("connection is not positioned at the start ",
                         "of the file, rewinding it"), immediate.=TRUE)
            seek(filexp, where=0)
        }
    }
    .Call(read_gff_pragmas, filexp)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sniffGFFVersion()
###

.get_version_from_pragmas <- function(pragmas)
{
    attrcol_fmt <- attr(pragmas, "attrcol_fmt")
    idx <- grep("^##gff-version", pragmas)
    if (length(idx) == 0L) {
        if (is.null(attrcol_fmt))
            stop(wmsg("'attr(pragmas, \"attrcol_fmt\")' is NULL"))
        return(attrcol_fmt)
    }
    version <- sub("^##gff-version", "", pragmas[idx])
    version <- unique(version)
    if (length(version) > 1L) {
        warning(wmsg("more than one GFF version specified in the file, ",
                     "returning the first one"))
        version <- version[[1L]]
    }
    version <- suppressWarnings(as.integer(version))
    if (is.na(version)) {
        warning(wmsg("unrecognized GFF version specified in the file"))
        if (is.null(attrcol_fmt))
            stop(wmsg("'attr(pragmas, \"attrcol_fmt\")' is NULL"))
        return(attrcol_fmt)
    }
    version
}

sniffGFFVersion <- function(filepath)
{
    pragmas <- readGFFPragmas(filepath)
    .get_version_from_pragmas(pragmas)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GFFcolnames()
###
### Return the 9 standard GFF columns as specified at:
###   http://www.sequenceontology.org/resources/gff3.html
###

GFFcolnames <- function(GFF1=FALSE)
{
    if (!isTRUEorFALSE(GFF1))
        stop(wmsg("'GFF1' must be TRUE or FALSE"))
    .Call(gff_colnames, GFF1)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGFF()
###

### Returns 0L, 1L, 2L, or 3L.
.normarg_version <- function(version=0)
{
    if (is.character(version)) {
        ## For compatibility with "import" method for GFFFile objects.
        if (length(version) == 0L)
            return(0L)
        if (isSingleString(version)) {
            IMPORT_STYLE_VERSIONS <- c("", "1", "2", "3")
            m <- match(version, IMPORT_STYLE_VERSIONS)
            if (is.na(m))
                stop(wmsg("when a single string, 'version' must ",
                          "be \"\", \"1\", \"2\", or \"3\""))
            version <- m - 1L
            return(version)
        }
    }
    if (isSingleNumber(version)) {
        if (!is.integer(version))
            version <- as.integer(version)
        if (version < 0L || version > 3L)
            stop(wmsg("'version' must be 0, 1, 2, or 3"))
        return(version)
    }
    stop(wmsg("'version' must be a single number"))
}

.prepare_colmap_and_tags <- function(columns=NULL, tags=NULL, attrcol_fmt=0L)
{
    ## Check 'columns'.
    if (!(is.null(columns) || is.character(columns)))
        stop(wmsg("'columns' must be NULL or a character vector"))

    GFF_colnames <- GFFcolnames(attrcol_fmt == 1L)

    ## Check 'tags'.
    if (!is.null(tags)) {
        if (!is.character(tags))
            stop(wmsg("'tags' must be NULL or character vector"))
        if (any(is.na(tags)) || anyDuplicated(tags))
            stop(wmsg("'tags' cannot contain NAs or duplicates"))
        if (attrcol_fmt == 1L) {
            ## Move any GFF colname found in 'tags' to 'columns'.
            columns <- union(columns, intersect(tags, GFF_colnames))
            tags <- setdiff(tags, GFF_colnames)
            if (length(tags) != 0L)
                warning(wmsg("trying to extract tags from a GFF1 file"))
        }
    }

    ## Prepare 'colmap'.
    if (is.null(columns)) {
        colmap <- seq_along(GFF_colnames)
        if (attrcol_fmt != 1L) {
            stopifnot(GFF_colnames[[length(GFF_colnames)]] == "attributes")
            ## We don't load the "attributes" column unless the user
            ## requested no tags (i.e. by setting 'tags' to character(0)).
            if (!(is.character(tags) && length(tags) == 0L))
                colmap[[length(GFF_colnames)]] <- NA_integer_
        }
    } else if (is.character(columns)) {
        if (!all(columns %in% GFF_colnames)) {
            in1string <- paste0(GFF_colnames, collapse=", ")
            stop(wmsg("'columns' must contain valid GFF columns. ",
                  "Valid GFF columns are: ", in1string))
        }
        if (anyDuplicated(columns))
            stop(wmsg("'columns' cannot contain duplicates"))
        colmap <- match(GFF_colnames, columns)
    } else {
        stop(wmsg("'columns' must be NULL or a character vector"))
    }

    list(colmap=colmap, tags=tags)
}

.normarg_filter <- function(filter, attrcol_fmt=0L)
{
    GFF_colnames <- GFFcolnames(attrcol_fmt == 1L)
    if (is.null(filter))
        return(NULL)
    if (!is.list(filter))
        stop(wmsg("'filter' must be NULL or a named list"))
    filter_names <- names(filter)
    if (is.null(filter_names))
        stop(wmsg("'filter' must have names"))
    if (attrcol_fmt == 1L) {
        valid_filter_names <- GFF_colnames
    } else {
        valid_filter_names <- head(GFF_colnames, n=-1L)
    }
    if (!all(filter_names %in% valid_filter_names)) {
        in1string <- paste0(valid_filter_names, collapse=", ")
        if (attrcol_fmt == 1L) {
            excluding_note <- ""
        } else {
            excluding_note <- "(excluding \"attributes\")"
        }
        stop(wmsg("The names on 'filter' must be valid GFF columns ",
                  excluding_note, ". ",
                  "Valid 'filter' names: ", in1string))
    }
    if (anyDuplicated(filter_names))
        stop(wmsg("the names on 'filter' must be unique"))
    unname(filter[valid_filter_names])
}

### 'df' must be a data-frame-like object (typically an ordinary data frame or
### a DataFrame object).
.is_multi_tag <- function(df, ntag, attrcol_fmt=0L)
{
    if (ntag == 0L || attrcol_fmt != 3L)
        return(logical(ntag))
    multi_tags <- c("Parent", "Alias", "Note",
                    "Dbxref", "Ontology_term")
    sapply(seq_len(ntag) + ncol(df) - ntag,
           function(j)
               colnames(df)[[j]] %in% multi_tags ||
                   any(grepl(",", df[[j]], fixed=TRUE)))
}

### 'df' must be a data-frame-like object (typically an ordinary data frame or
### a DataFrame object). 'decode_idx' must be a non-empty integer vector
### indicating which columns to decode. The columns to decode must be character
### vectors.
.url_decode_cols <- function(df, decode_idx)
{
    decoded_cols <- lapply(setNames(decode_idx, colnames(df)[decode_idx]),
                           function(j)
                             urlDecode(df[[j]], na.strings=NA_character_))
    df[decode_idx] <- decoded_cols
    df
}

### 'df' must be a data-frame-like object (typically an ordinary data frame or
### a DataFrame object). 'split_idx' must be a non-empty integer vector
### indicating which columns to split. The columns to split must be character
### vectors. Split values are passed thru urlDecode() unless 'raw_data' is
### TRUE. Always returns a DataFrame.
.strsplit_cols <- function(df, split_idx, raw_data)
{
    split_cols <- lapply(setNames(split_idx, colnames(df)[split_idx]),
        function(j) {
            col <- df[[j]]
            ## Probably the most efficient way to create an empty CharacterList
            ## of arbitrary length.
            split_col <- relist(character(0),
                                PartitioningByEnd(rep.int(0L, length(col))))
            not_na <- !is.na(col)
            tmp <- strsplit(col[not_na], ",", fixed=TRUE)
            split_col[not_na] <- CharacterList(tmp)
            if (raw_data)
                return(split_col)
            relist(urlDecode(unlist(split_col)), split_col)
        })
    ## Surprisingly sticking the CharacterList cols back into 'df' works
    ## even if 'df' is an ordinary data frame!
    df[split_idx] <- split_cols
    ans <- DataFrame(df, check.names=FALSE)
    ## "show" method for DataFrame is broken if some colnames are the empty
    ## string so we rename this column (in our case, we know there can only
    ## be one).
    m <- match("", colnames(ans))
    if (!is.na(m))
        colnames(ans)[m] <- "__empty_tag__"
    ans
}

readGFF <- function(filepath, version=0, columns=NULL, tags=NULL,
                    filter=NULL, nrows=-1, raw_data=FALSE)
{
    ## Check 'filepath'.
    filexp <- .make_filexp_from_filepath(filepath)
    if (is(filexp, "connection")) {
        if (!isOpen(filexp)) {
            open(filexp)
            on.exit(close(filexp))
        }
        if (seek(filexp) != 0) {
            warning(wmsg("connection is not positioned at the start ",
                         "of the file, rewinding it"), immediate.=TRUE)
            seek(filexp, where=0)
        }
    }

    ## Check 'version'.
    version <- .normarg_version(version)

    ## Get pragmas lines.
    pragmas <- .Call(read_gff_pragmas, filexp)

    ## Rewind file.
    if (is(filexp, "connection")) {
        seek(filexp, where=0)
    } else {
        XVector:::rewind_filexp(filexp)
    }

    if (version == 0L) {
        attrcol_fmt <- .get_version_from_pragmas(pragmas)
    } else {
        attrcol_fmt <- version
    }

    ## Prepare 'colmap' and normalize 'tags'.
    colmap_and_tags <- .prepare_colmap_and_tags(columns, tags, attrcol_fmt)
    colmap <- colmap_and_tags$colmap
    tags <- colmap_and_tags$tags

    ## Normalize 'filter'.
    filter <- .normarg_filter(filter, attrcol_fmt)

    ## Normalize 'nrows'.
    if (!isSingleNumber(nrows))
        stop(wmsg("'nrows' must be a single number"))
    if (!is.integer(nrows))
        nrows <- as.integer(nrows)

    ## Check 'raw_data'.
    if (!isTRUEorFALSE(raw_data))
        stop(wmsg("'raw_data' must be TRUE or FALSE"))

    ## 1st pass.
    scan_ans <- .Call(scan_gff, filexp, attrcol_fmt, tags, filter, nrows)
    if (is.null(tags))
        tags <- scan_ans[[1L]]
    nrows <- scan_ans[[2L]]

    ## Rewind file.
    if (is(filexp, "connection")) {
        seek(filexp, where=0)
    } else {
        XVector:::rewind_filexp(filexp)
    }

    ## 2nd pass: return 'ans' as an ordinary data frame.
    ans <- .Call(load_gff, filexp, attrcol_fmt, tags, filter,
                           nrows, pragmas,
                           colmap, raw_data)
    ncol0 <- attr(ans, "ncol0")
    ntag <- attr(ans, "ntag")          # should be the same as 'length(tags)'

    ## Post-process standard GFF cols.
    if (!raw_data) {
        if (attrcol_fmt == 1L) {
            #factor_colnames <- c("seqid", "source", "type", "strand", "group")
            factor_colnames <- c("seqid", "source", "type", "group")
        } else {
            #factor_colnames <- c("seqid", "source", "type", "strand")
            factor_colnames <- c("seqid", "source", "type")
        }
        m <- match(factor_colnames, head(colnames(ans), n=ncol0))
        m <- m[!is.na(m)]
        factor_cols <- lapply(setNames(m, colnames(ans)[m]),
                              function(j)
                                factor(ans[[j]], levels=unique(ans[[j]])))
        ans[m] <- factor_cols
    }

    ## Post-process tags.
    if (ntag != 0L) {
        is_multi_tag <- .is_multi_tag(ans, ntag, attrcol_fmt)
        if (!raw_data) {
            decode_idx <- which(!is_multi_tag) + ncol0
            if (length(decode_idx) != 0L)
                ans <- .url_decode_cols(ans, decode_idx)
        }
        split_idx <- which(is_multi_tag) + ncol0
        if (length(split_idx) != 0L) {
            ## Returns 'ans' as a DataFrame.
            ans <- .strsplit_cols(ans, split_idx, raw_data)
        }
    }

    ## 'ans' could have lost its readGFF-specific attributes (e.g. if it was
    ## turned into a DataFrame), so we restore them and cross our fingers that
    ## they won't clash with the DataFrame slots the day the internals of
    ## DataFrame objects happen to change (very unlikely though).
    if (is.null(attr(ans, "pragmas")))
        attr(ans, "pragmas") <- pragmas
    if (is.null(attr(ans, "attrcol_fmt")))
        attr(ans, "attrcol_fmt") <- attrcol_fmt
    if (is.null(attr(ans, "ncol0")))
        attr(ans, "ncol0") <- ncol0
    if (is.null(attr(ans, "ntag")))
        attr(ans, "ntag") <- ntag
    if (is.null(attr(ans, "raw_data")))
        attr(ans, "raw_data") <- raw_data
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readGFFAsGRanges()
###
### Used in "import" method for GFFFile objects.
### Not exported (user should use import()).
###

### sequence-region => Seqinfo -- by Michael
.parseSequenceRegionsAsSeqinfo <- function(lines) {
    sr <- grep("##sequence-region", lines, value=TRUE)
    srcon <- file()
    on.exit(close(srcon))
    writeLines(sr, srcon)
    srt <- read.table(srcon, comment.char="",
                      colClasses=list(NULL, "character", "integer",
                          "integer"))
    if (any(srt[[2L]] != 1L)) {
        warning("One or more ##sequence-region directives do not start at 1. ",
                "The assumptions made by 'sequenceRegionsAsSeqinfo=TRUE' ",
                "have been violated.")
    }
    Seqinfo(srt[[1L]], srt[[3L]])
}

### -- by Michael
.parseSpeciesAsMetadata <- function(lines) {
    species <- unique(grep("##species", lines, fixed=TRUE, value=TRUE))
    if (length(species) > 1L) {
        stop("multiple species definitions found")
    }
    metadata <- list()
    if (length(species) == 1L) {
        species <- sub("##species ", "", species, fixed=TRUE)
        if (isNCBISpeciesURL(species)) {
            ncbiError <- function(e) {
                warning("failed to retrieve organism information from NCBI")
            }
            metadata <- tryCatch(metadataFromNCBI(species), error = ncbiError)
        }
    }
    metadata
}

readGFFAsGRanges <- function(filepath, version=0, colnames=NULL, filter=NULL,
                             genome=NA,
                             sequenceRegionsAsSeqinfo=FALSE,
                             speciesAsMetadata=FALSE)
{
    if (!isSingleStringOrNA(genome))
        stop(wmsg("'genome' must be a single string or NA"))
    if (!isTRUEorFALSE(sequenceRegionsAsSeqinfo))
        stop(wmsg("'sequenceRegionsAsSeqinfo' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(speciesAsMetadata))
        stop(wmsg("'speciesAsMetadata' must be TRUE or FALSE"))

    ## Read as data frame.
    if (is.null(colnames)) {
        df <- readGFF(filepath, version=version, filter=filter)
    } else {
        if (!is.character(colnames))
            stop(wmsg("'colnames' must be a character vector"))
        ## Split 'colnames' between 'columns' and 'tags'.
        GFF_colnames <- GFFcolnames()
        columns <- intersect(colnames, GFF_colnames)
        tags <- setdiff(colnames, GFF_colnames)
        core_columns <- c("seqid", "start", "end", "strand")
        columns <- union(columns, core_columns)
        df <- readGFF(filepath, version=version,
                      columns=columns, tags=tags, filter=filter)
    }

    ## Get 'ans_seqinfo' from pragmas.
    pragmas <- attr(df, "pragmas")
    attrcol_fmt <- attr(df, "attrcol_fmt")
    if (sequenceRegionsAsSeqinfo && attrcol_fmt == 3L) {
        ans_seqinfo <- .parseSequenceRegionsAsSeqinfo(pragmas)
    } else if (is.na(genome)) {
        ans_seqinfo <- NULL
    } else {
        ans_seqinfo <- seqinfoForGenome(genome)
    }

    ## Get 'ans_metadata' from pragmas.
    if (speciesAsMetadata) {
        ans_metadata <- .parseSpeciesAsMetadata(pragmas)
    } else {
        ans_metadata <- list()
    }

    ## Turn data frame into GRanges.
    ## TODO: Maybe we should be able to pass the metadata to
    ## makeGRangesFromDataFrame()?
    if (is.null(colnames)) {
        ans <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE,
                                            seqinfo=ans_seqinfo,
                                            seqnames.field="seqid")
    } else {
        ans <- makeGRangesFromDataFrame(df, seqinfo=ans_seqinfo,
                                            seqnames.field="seqid")
        mcols(ans) <- df[ , colnames, drop=FALSE]
    }
    genome(ans) <- genome
    metadata(ans) <- ans_metadata
    ans
}

