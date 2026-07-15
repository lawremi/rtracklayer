### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Refer methods
###
is_valid := new_generic("ref")

method(is_valid, Refer) <- function(ref) {
    child_value <- prop(ref@child, ref@child_key)
    if (is.null(child_value))
        return(NULL)
    missing <- setdiff(child_value, ref@valid)
    if (length(missing) > 0L) {
        ctx <- if (is.null(ref@context)) ref@child_key
               else paste0(ref@child_key, " (", ref@context, ")")
        paste0("invalid ", ctx, ": ",
               paste0("'", missing, "'", collapse = ", "),
               " not in {",
               paste0(ref@valid, collapse = ", "), "}")
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Class dispatch tables
###

## key presence → class mapping for pick_class()
genome_choices <- list("twoBitPath" = Assembly,
                       "!_default" = Genome)

track_choices <- list("superTrack" = SuperTrack,
                      "compositeTrack" = CompositeTrack,
                      "container" = OverlayTrack,
                      "view" = View,
                      "!_default" = Track)

## type value → subclass mapping for refine_track_class()
type_classes <- list("bam" = BamTrack,
                     "bigBarChart" = BigBarChartTrack,
                     "bigBed" = BigBedTrack,
                     "bigChain" = BigChainTrack,
                     "bigGenePred" = BigGenePredTrack,
                     "bigInteract" = BigInteractTrack,
                     "bigLolly" = BigLollyTrack,
                     "bigMaf" = BigMafTrack,
                     "bigNarrowPeak" = BigNarrowPeakTrack,
                     "bigPsl" = BigPslTrack,
                     "bigWig" = BigWigTrack,
                     "halSnake" = HalSnakeTrack,
                     "hic" = HicTrack,
                     "vcfTabix" = VcfTabixTrack,
                     "vcfPhasedTrio" = VcfPhasedTrioTrack)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Reading hub files
###

read_hub := new_generic(c("uri", "filename"))

## Build a list of track objects from a data frame of track stanzas.
## Shared by both the useOneFile and multi-file branches of read_hub.
build_tracks <- function(track_stanzas) {
    unlist(lapply(track_stanzas, function(track) {
        cls <- pick_class(track, track_choices)
        cls <- refine_track_class(track, cls)
        stanzas_to(track, cls)
    }), recursive = FALSE)
}

## Load per-genome metadata (tagStorm first, then tab-separated) if the
## referenced file actually exists under `uri`. Returns NULL when neither
## file is present, which is the common case for minimal hubs.
load_metadata <- function(uri, genome) {
    metadb <- genome@metaDb
    if (!is.null(metadb) && file.exists(append_to_uri(uri, metadb)))
        return(parse_metadb(uri, metadb))
    metatab <- genome@metaTab
    if (!is.null(metatab) && file.exists(append_to_uri(uri, metatab)))
        return(parse_metatab(uri, metatab))
    NULL
}

## Load groups.txt for a genome if present. Returns a named list of Group
## objects, or NULL. Only Assembly carries a `groups` slot; plain Genome
## does not, so this is a no-op for non-Assembly genomes.
load_groups <- function(uri, genome) {
    if (!S7_inherits(genome, Assembly)) return(NULL)
    groups_file <- genome@groups
    if (is.null(groups_file)) return(NULL)
    path <- append_to_uri(uri, groups_file)
    if (!file.exists(path)) return(NULL)
    group_df <- parse_file(uri, groups_file, stanza_keys = STANZA_KEYS_GROUPS)
    stanzas_to(split(group_df, group_df$group), Group)
}

method(read_hub, list(class_character, class_character)) <-
    function(uri, filename) {
    stanza_df <- parse_file(uri, filename)
    stanza_list <- split(stanza_df, stanza_df$group)

    ## Single-stanza hub.txt → multi-file layout. Follow genomesFile to
    ## fetch genome stanzas, then for each genome fetch its trackDb plus
    ## optional groups/metadata files.
    if (length(stanza_list) == 1L) {
        hub <- stanzas_to(stanza_list[[1L]], Hub)[[1L]]
        genomes_df <- parse_file(uri, hub@genomesFile)
        genome_stanzas <- split(genomes_df, genomes_df$group)
        genomes <- stanzas_to(
            genome_stanzas,
            pick_class(genome_stanzas[[1L]], genome_choices)
        )
        tracks <- unlist(lapply(genomes, function(genome) {
            track_df <- parse_file(uri, genome@trackDb)
            build_tracks(split(track_df, track_df$group))
        }), recursive = FALSE)
        ## Metadata and groups come from the first genome that has them.
        ## Multi-genome hubs that carry distinct metadata per genome are
        ## a future extension.
        metadata <- NULL
        groups <- NULL
        for (genome in genomes) {
            if (is.null(metadata)) metadata <- load_metadata(uri, genome)
            if (is.null(groups)) groups <- load_groups(uri, genome)
        }
        return(TrackHub(hub = hub, genomes = genomes, tracks = tracks,
                        metadata = metadata, groups = groups))
    }

    ## useOneFile layout: hub stanza, then genome stanza(s), then tracks.
    hub_stanza <- stanza_list[[1L]]
    genome_stanza <- stanza_list[[2L]]
    track_stanzas <- stanza_list[-c(1L, 2L)]
    TrackHub(
        hub = stanzas_to(hub_stanza, Hub)[[1L]],
        genomes = stanzas_to(
            genome_stanza,
            pick_class(genome_stanza, genome_choices)
        ),
        tracks = build_tracks(track_stanzas)
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Stanza parsing (hub.txt, genomes.txt, trackDb.txt, groups.txt)
###

append_to_uri <- function(uri, filename) {
    if (uriIsWritable(uri))
        file.path(uri, filename)
    else
        paste0(gsub(RE_TRAILING_SLASH, "/", uri), filename)
}

join_continuations <- function(lines) {
    if (length(lines) == 0L) return(lines)
    block_starts <- c(TRUE, !grepl(RE_CONTINUATION, lines[-length(lines)]))
    groups <- cumsum(block_starts)
    cleaned <- sub(RE_CONTINUATION, "", lines)
    is_cont <- !block_starts
    cleaned[is_cont] <- trimws(cleaned[is_cont])
    unname(tapply(cleaned, groups, paste, collapse = " "))
}

## Default stanza-boundary keys for the main hub files. Adjust via the
## stanza_keys argument to parse_file() when reading files whose stanzas
## start with a different key (e.g. groups.txt starts each stanza with
## `name`).
STANZA_KEYS_HUB <- c("track", "genome")
STANZA_KEYS_GROUPS <- c("name")

parse_file <- function(uri, filename, stanza_keys = STANZA_KEYS_HUB) {
    file_path <- append_to_uri(uri, filename)
    lines <- tryCatch(
        readLines(file_path, warn = FALSE),
        error = function(e) {
            stop("failed to read '", file_path, "': ",
                 conditionMessage(e), call. = FALSE)
        }
    )
    lines <- join_continuations(lines)
    lines <- trimws(lines)
    lines <- lines[lines != "" & !grepl(RE_COMMENT, lines)]
    if (length(lines) == 0L)
        stop("no content found in '", file_path, "'", call. = FALSE)
    lines <- gsub(RE_TAB, " ", lines)
    key <- sub(RE_KEY_EXTRACT, "\\1", lines)
    value <- sub(RE_VALUE_EXTRACT, "\\1", lines)
    value[!grepl(RE_HAS_SPACE, lines)] <- ""
    stanza_pattern <- paste0("^(", paste(stanza_keys, collapse = "|"), ")\\b")
    data.frame(
        group = cumsum(grepl(stanza_pattern, lines)),
        key = key,
        value = value
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Metadata parsing (tagStorm & tab-separated)
###

parse_metatab <- function(uri, filename) {
    file_path <- append_to_uri(uri, filename)
    MetaData(table = read.table(file_path, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE))
}

parse_metadb <- function(uri, filename) {
    file_path <- append_to_uri(uri, filename)
    lines <- tryCatch(
        readLines(file_path, warn = FALSE),
        error = function(e) {
            stop("failed to read '", file_path, "': ",
                 conditionMessage(e), call. = FALSE)
        }
    )
    MetaData(table = parse_tagstorm(lines))
}

parse_tagstorm <- function(lines) {
    lines <- lines[!grepl(RE_COMMENT_WS, lines)]
    if (length(lines) == 0L)
        return(data.frame())

    raw_indent <- nchar(sub(RE_LEADING_SPACES, "\\1", lines))
    is_blank <- grepl(RE_BLANK, lines)
    trimmed <- trimws(lines)
    key <- sub(RE_KEY_EXTRACT_OPT, "\\1", trimmed)
    value <- sub(RE_VALUE_EXTRACT_OPT, "", trimmed)

    ## Stack of key-value pairs keyed by indent depth.
    ## Blank lines signal sibling boundaries — prune current depth and deeper.
    ## Consecutive lines at the same depth accumulate into the same stanza.
    stack <- list()
    rows <- list()

    for (i in seq_along(lines)) {
        if (is_blank[i]) next
        depth <- raw_indent[i]
        depth_key <- as.character(depth)

        preceded_by_blank <- i > 1L && is_blank[i - 1L]
        if (preceded_by_blank)
            stack <- stack[as.integer(names(stack)) < depth]
        else
            stack <- stack[as.integer(names(stack)) <= depth]

        if (is.null(stack[[depth_key]]))
            stack[[depth_key]] <- character()
        stack[[depth_key]][[key[i]]] <- value[i]

        ## "meta" key marks a complete leaf — flatten with inheritance
        if (key[i] == "meta") {
            flat <- character()
            for (d in as.character(sort(as.integer(names(stack)))))
                flat[names(stack[[d]])] <- stack[[d]]
            rows <- c(rows, list(flat))
        }
    }

    if (length(rows) == 0L)
        return(data.frame())

    all_keys <- unique(unlist(lapply(rows, names)))
    df <- do.call(rbind, lapply(rows, function(row) {
        vals <- setNames(rep(NA_character_, length(all_keys)), all_keys)
        vals[names(row)] <- row
        as.data.frame(as.list(vals), stringsAsFactors = FALSE)
    }))
    rownames(df) <- df$meta
    df
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Stanza → S7 object coercion
###

pick_class <- function(stanza, choices) {
    for (key in names(choices)) {
        should_pick <- !is.null(stanza[[key]])
        if (grepl(RE_NEGATION_PREFIX, key))
            should_pick <- is.null(stanza[[key]])
        if (should_pick)
            return(choices[[key]])
    }
    keys_present <- if (is.data.frame(stanza)) stanza$key else names(stanza)
    stop("unable to determine class for stanza with keys: ",
         paste(keys_present, collapse = ", "), call. = FALSE)
}

refine_track_class <- function(stanza, class) {
    if (!identical(class, Track)) return(class)
    type_val <- if (is.data.frame(stanza)) {
        vals <- stanza$value[stanza$key == "type"]
        if (length(vals) == 0L) return(class)
        vals[1L]
    } else stanza[["type"]]
    if (is.null(type_val)) return(class)
    base <- sub(RE_TYPE_BASE, "", trimws(type_val))
    type_classes[[base]] %||% class
}

get_prop_class <- function(prop) {
    cls <- prop$class
    if (inherits(cls, "S7_union")) {
        classes <- Filter(Negate(is.null), cls$classes)
        base_classes <- vapply(classes, `[[`, character(1L), "class")
        if (all(c("integer", "double") %in% base_classes))
            return("numeric")
        return(classes[[1L]]$class)
    }
    cls$class
}

## Routing rules for dotted dynamic settings.
## Each entry: target slot name → vector of accepted dotted prefixes.
DOTTED_PREFIX_ROUTES <- list(
    filters = c("filter.", "filterText.", "filterValues.", "filterLabel."),
    highlights = c("highlight.", "highlightText.", "highlightValues."),
    decorators = c("decorator."),
    yAxisLabels = c("yAxisLabel.")
)

route_dotted_args <- function(args, valid_args) {
    routed <- list()
    consumed <- logical(length(args))
    for (slot in names(DOTTED_PREFIX_ROUTES)) {
        if (!(slot %in% valid_args)) next
        prefixes <- DOTTED_PREFIX_ROUTES[[slot]]
        pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
        matched <- grepl(pattern, names(args))
        if (any(matched)) {
            routed[[slot]] <- as.list(args[matched])
            consumed <- consumed | matched
        }
    }
    list(routed = routed, remaining = args[!consumed])
}

extract_args <- function(stanza, valid_args) {
    args <- setNames(stanza$value, stanza$key)
    routed <- route_dotted_args(args, valid_args)
    direct <- routed$remaining[intersect(names(routed$remaining), valid_args)]
    c(as.list(direct), routed$routed)
}

coerce_args <- function(args, props_class) {
    result <- list()
    for (prop in names(args)) {
        cls <- props_class[prop]
        value <- args[[prop]]
        if (!is(value, cls) && !is(value, "list")) {
            as_fun <- get(paste0("as.", cls), mode = "function",
                          inherits = TRUE)
            value <- as_fun(value)
        }
        result[[prop]] <- value
    }
    result
}

## Warn (but don't error) when a label property exceeds its UCSC-spec
## character cap. Called once per stanza, post-coercion, by stanzas_to().
check_label_caps <- function(args, stanza_id) {
    for (name in intersect(names(args), names(LABEL_CAPS))) {
        cap <- LABEL_CAPS[[name]]
        value <- args[[name]]
        if (is.null(value) || nchar(value) <= cap) next
        warning(stanza_id, ": ", name, " exceeds ", cap,
                " character limit (", nchar(value), " chars): '",
                value, "'", call. = FALSE)
    }
}

construct_objects <- function(coerced_args, class, id) {
    objects <- lapply(coerced_args, function(args) do.call(class, args))
    names(objects) <- vapply(objects, function(obj) prop(obj, id),
                             character(1L))
    objects
}

stanzas_to <- function(stanza_list, class, id = NULL) {
    if (is.data.frame(stanza_list)) stanza_list <- list(stanza_list)
    valid_args <- names(class@properties)
    if (is.null(id)) id <- valid_args[1L]
    props_class <- vapply(class@properties, get_prop_class, character(1L))
    args_list <- lapply(stanza_list, extract_args, valid_args)
    coerced <- lapply(args_list, coerce_args, props_class)
    for (args in coerced)
        check_label_caps(args, stanza_id = args[[id]] %||% "<unnamed>")
    construct_objects(coerced, class, id)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show methods
###

strip_prefix <- function(x) sub(RE_STRIP_NUM_PREFIX, "", x)

track_class_label <- function(track) {
    cls <- S7_class(track)
    if (is.null(cls)) return("Track")
    cls@name
}

method(print, TrackHub) <- function(x, ...) {
    hub <- x@hub
    header <- sprintf("TrackHub: %s", hub@shortLabel)
    email_line <- sprintf("  email: %s", hub@email)
    onefile <- if (!is.null(hub@useOneFile) && hub@useOneFile == "on")
        "  useOneFile: on"

    genome_names <- vapply(x@genomes, function(g) g@genome, character(1L))
    genome_line <- sprintf("  genomes(%d): %s",
                           length(genome_names),
                           paste(genome_names, collapse = ", "))

    track_names <- strip_prefix(names(x@tracks))
    track_types <- vapply(x@tracks, function(t) track_class_label(t),
                          character(1L))
    track_header <- sprintf("  tracks(%d):", length(x@tracks))
    track_lines <- sprintf("    - %s [%s]", track_names, track_types)

    meta_line <- if (length(x@metadata) > 0L)
        sprintf("  metadata(%d)", length(x@metadata))
    group_line <- if (length(x@groups) > 0L)
        sprintf("  groups(%d): %s", length(x@groups),
                paste(vapply(x@groups, function(g) g@name, character(1L)),
                      collapse = ", "))

    cat(c(header, email_line, onefile, genome_line,
          track_header, track_lines, meta_line, group_line),
        sep = "\n")
    invisible(x)
}
