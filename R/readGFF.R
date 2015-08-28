### =========================================================================
### readGFF()
### -------------------------------------------------------------------------


### Return the 9 standard GFF columns as specified at:
###   http://www.sequenceontology.org/resources/gff3.html
GFFcolnames <- function() .Call(gff_colnames)

### 'df' must be a data-frame-like object (typically an ordinary data frame or
### a DataFrame object). 'decode_idx' must be a non-empty integer vector
### indicating which columns to decode. The columns to decode must be character
### vectors.
.urlDecodeCols <- function(df, decode_idx)
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
.strsplitCols <- function(df, split_idx, raw_data)
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

### Does NOT work on a connection object.
readGFF <- function(filepath, columns=NULL, tags=NULL,
                    filter=NULL, raw_data=FALSE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string"))
    if (!isTRUEorFALSE(raw_data))
        stop(wmsg("'raw_data' must be TRUE or FALSE"))
    filexp <- XVector:::open_input_files(filepath)[[1L]]

    ## Check 'tags'.
    if (!is.null(tags)) {
        if (!is.character(tags))
            stop(wmsg("'tags' must be NULL or character vector"))
        if (any(is.na(tags)) || anyDuplicated(tags))
            stop(wmsg("'tags' cannot contain NAs or duplicates"))
    }

    ## Prepare 'colmap'.
    GFF_colnames <- GFFcolnames()
    stopifnot(GFF_colnames[[length(GFF_colnames)]] == "attributes")
    if (is.null(columns)) {
        colmap <- seq_along(GFF_colnames)
        ## We don't load the "attributes" column unless the user requested no
        ## tags (i.e. by setting 'tags' to character(0)).
        if (!(is.character(tags) && length(tags) == 0L))
            colmap[[length(GFF_colnames)]] <- NA_integer_
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
        stop(wmsg("'columns' must be NULL or character vector"))
    }

    ## Normalize 'filter'.
    if (!is.null(filter)) {
        if (!is.list(filter))
            stop("'filter' must be NULL or a named list")
        filter_names <- names(filter)
        if (is.null(filter_names))
            stop("'filter' must have names")
        valid_filter_names <- head(GFF_colnames, n=-1L)
        if (!all(filter_names %in% valid_filter_names)) {
            in1string <- paste0(valid_filter_names, collapse=", ")
            stop(wmsg("The names on 'filter' must be valid GFF columns ",
                      "(excluding \"attributes\"). ",
                      "Valid names on 'filter': ", in1string))
        }
        if (anyDuplicated(filter_names))
            stop(wmsg("names on 'filter' must be unique"))
        filter <- filter[valid_filter_names]
    }

    ## 1st pass.
    scan_ans <- .Call(scan_gff, filexp, tags, filter)
    if (is.null(tags))
        tags <- scan_ans[[1L]]
    ans_nrow <- scan_ans[[2L]]
    attrcol_fmt <- scan_ans[[3L]]
    pragmas <- scan_ans[[4L]]

    ## 2nd pass: return 'ans' as an ordinary data frame.
    ans <- .Call(load_gff, filexp, tags, filter,
                           ans_nrow, attrcol_fmt, pragmas,
                           colmap, raw_data)
    ncol0 <- attr(ans, "ncol0")
    ntag <- attr(ans, "ntag")          # should be the same as 'length(tags)'
    raw_data <- attr(ans, "raw_data")  # should be the same as user-supplied
    pragmas <- attr(ans, "pragmas")

    ## Post-process standard GFF cols.
    #factor_colnames <- c("seqid", "source", "type", "strand")
    factor_colnames <- c("seqid", "source", "type")
    m <- match(factor_colnames, head(colnames(ans), n=ncol0))
    m <- m[!is.na(m)]
    factor_cols <- lapply(setNames(m, colnames(ans)[m]),
                          function(j)
                            factor(ans[[j]], levels=unique(ans[[j]])))
    ans[m] <- factor_cols

    ## Post-process tags.
    if (ntag != 0L) {
        multi_tags <- c("Parent", "Alias", "Note", "Dbxref", "Ontology_term")
        is_multi_tag <- sapply(seq_len(ntag) + ncol0,
                               function(j)
                                 colnames(ans)[[j]] %in% multi_tags ||
                                   any(grepl(",", ans[[j]], fixed=TRUE)))
        if (!raw_data) {
            decode_idx <- which(!is_multi_tag) + ncol0
            if (length(decode_idx) != 0L)
                ans <- .urlDecodeCols(ans, decode_idx)
        }
        split_idx <- which(is_multi_tag) + ncol0
        if (length(split_idx) != 0L) {
            ## Returns 'ans' as a DataFrame.
            ans <- .strsplitCols(ans, split_idx, raw_data)
        }
    }

    ## 'ans' could have lost its readGFF-specific attributes (e.g. if it was
    ## turned into a DataFrame), so we restore them and cross our fingers that
    ## they won't clash with the DataFrame slots the day the internals of
    ## DataFrame objects happen to change (very unlikely though).
    if (is.null(attr(ans, "ncol0")))
        attr(ans, "ncol0") <- ncol0
    if (is.null(attr(ans, "ntag")))
        attr(ans, "ntag") <- ntag
    if (is.null(attr(ans, "raw_data")))
        attr(ans, "raw_data") <- raw_data
    if (is.null(attr(ans, "pragmas")))
        attr(ans, "pragmas") <- pragmas
    ans
}

