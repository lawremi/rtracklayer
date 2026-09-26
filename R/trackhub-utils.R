## taken from lawremi/wizrd
`:=` <- function(x, y) {
    sym <- substitute(x)
    stopifnot(is.name(sym))
    call <- substitute(y)
    stopifnot(is.call(call))

    nm <- deparse(sym)
    call$name <- nm

    assign(nm, eval(call, parent.frame()), parent.frame())
}

nullable <- function(prop) {
    if (is.null(prop) || S7:::is_foundation_class(prop))
        prop <- new_property(prop)
    prop$class <- new_union(NULL, prop$class)
    if (identical(prop$default, missing_name()))
        prop["default"] <- list(NULL)
    prop
}

assert_scalar <- function(scalar, class, arg = deparse(substitute(scalar)))
{
    if (length(scalar) != 1L || !S7:::class_inherits(scalar, class)) {
        type_name <- if (identical(class, class_numeric)) {
            "numeric"
        } else if (inherits(class, "S7_union"))
            paste(class$classes, collapse = " | ")
        else class$class
        msg <- sprintf("`%s` must be a single %s value", arg, type_name)
        stop(msg, call. = FALSE)
    }
    if (is.na(scalar)) {
        msg <- sprintf("`%s` must not be NA", arg)
        stop(msg, call. = FALSE)
    }
}

missing_name <- function() alist(x=)[[1L]]

scalar <- function(x, ..., validator = x$validator, default = x$default,
                   choices = NULL)
{
    if (S7:::is_foundation_class(x))
        x <- new_property(x, ...)
    stopifnot(inherits(x, "S7_property"))
    x$default <- if (is.null(default)) {
        if (inherits(x$class, "S7_union") && is.null(x$class$classes[[1L]]))
            NULL
        else if (length(choices) > 0L)
            choices[1L]
        else missing_name()
    } else {
        if (!is.language(default))
            assert_scalar(default, x$class)
        default
    }
    force(validator)
    x$validator <- function(value) {
        if (is.null(value))
            return(NULL)
        c(if (length(value) != 1L || is.na(value))
            "must be of length one and not missing",
          if (!is.null(choices) && !all(value %in% choices))
                paste("contains values not in", deparse(choices)),
          if (!is.null(validator))
              validator(value)
          )
    }
    class(x) <- c("scalar_S7_property", class(x))
    x
}

class_object <- function(x) {
    if (is.null(x))
        return(NULL)
    S7_class(x) %||% as_class(getClassDef(class(x)[1L])) %||%
        new_S3_class(class(x))
}

zero_row_data_frame <- function(col.names) {
    data.frame(matrix(nrow = 0L, ncol = length(col.names))) |>
        setNames(col.names)
}

new_data_frame_property <- function(..., validator = NULL,
                                    col.names = colnames(prototype),
                                    default = substitute(prototype) %||%
                                        zero_row_data_frame(col.names),
                                    prototype = NULL)
{
    types <- lapply(prototype, class_object)
    prop <- new_property(class_data.frame, ..., validator = function(value) {
        c(if (!is.null(col.names) &&
                  !identical(colnames(value), col.names))
            paste("colnames() must be", deparse(col.names))
          else if (!is.null(prototype)) {
                wrong_type <- !mapply(inherits, value, types)
                if (any(wrong_type))
                    paste(colnames(value)[wrong_type], "must be a",
                          vapply(types[wrong_type], S7:::class_desc,
                                 character(1L)),
                          collapse = ", ")
          },
        if (!is.null(validator))
            validator(value)
        )
    }, default = default)
    prop$prototype <- prototype
    prop$col.names <- col.names
    class(prop) <- c("data_frame_S7_property", class(prop))
    prop
}

list_of <- function(class, ...) {
    new_list_property(of = class, ...)
}

new_list_property <- function(..., validator = NULL,
                              default = if (isTRUE(named))
                                  quote(setNames(list(), character()))
                              else quote(list()),
                              of = class_any, named = NA,
                              min_length = 0L, max_length = Inf) {
    prop <- new_property(class_list, ..., validator = function(value) {
        if (is.null(value))
            return(NULL)
        c(
            if (!identical(of, class_any) &&
                    !all(vapply(value, S7:::class_inherits, logical(1L), of)))
                paste("must only contain elements of class",
                      S7:::class_desc(of)),
            if (!is.null(of_validator)) {
                msgs <- unlist(lapply(value, of_validator))
                if (length(msgs) > 0L) {
                    paste("element(s) failed validation:",
                          paste0("'", unique(msgs), "'", collapse = ", "))
                }
            },
            if (isTRUE(named) && is.null(names(value)))
                "must have names",
            if (identical(named, FALSE) && !is.null(names(value)))
                "must not have names",
            if (length(value) < min_length || length(value) > max_length)
                paste0("must have length in [", min_length, ", ", max_length,
                       "]"),
            if (!is.null(validator))
                validator(value)
        )
    }, default = default)
    prop$of <- of
    if (inherits(of, "S7_property")) {
        of_validator <- of$validator
        of <- of$class
    } else {
        of_validator <- NULL
    }
    prop$named <- named
    class(prop) <- c("list_S7_property", class(prop))
    prop
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constants
###

HUB_FILE_NAME <- "hub.txt"
GENOME_FILE_NAME <- "genomes.txt"
METADB_FILE_NAME <- "tagStormFile.txt"
METATAB_FILE_NAME <- "tabSeparatedFile.txt"
GROUP_FILE_NAME <- "groups.txt"
TRACKDB_FILE_NAME <- "trackDb.txt"
TRACK_TYPES <- c("bam", "bigBed", "bigBarChart", "bigChain", "bigGenePred",
                 "bigInteract", "bigLolly", "bigNarrowPeak", "bigMaf",
                 "bigPsl", "bigWig", "hic", "halSnake", "vcfTabix",
                 "vcfPhasedTrio")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Regex
###

## URI / path helpers
RE_TRAILING_SLASH    <- "/*$"
RE_VALID_SCHEME      <- "^[A-Za-z][A-Za-z0-9+.-]*://"
RE_HTTP_FTP_SCHEME   <- "^(https?|ftp)://"

## Line-level parsing
RE_CONTINUATION      <- "\\s*\\\\\\s*$"
RE_COMMENT           <- "^#"
RE_COMMENT_WS        <- "^\\s*#"
RE_BLANK             <- "^\\s*$"
RE_HAS_SPACE         <- "\\s"
RE_WHITESPACE        <- "\\s+"
RE_LEADING_SPACES    <- "^( *).*"
RE_KEY_EXTRACT       <- "^(\\S+)\\s.*$"
RE_KEY_EXTRACT_OPT   <- "^(\\S+)\\s?.*$"
RE_VALUE_EXTRACT     <- "^\\S+\\s(.*)$"
RE_VALUE_EXTRACT_OPT <- "^\\S+\\s?"
RE_TAB               <- "\t"
RE_NEGATION_PREFIX   <- "^!"
RE_TYPE_BASE         <- "\\s.*$"
RE_STRIP_NUM_PREFIX  <- "^[0-9]+\\."

## Property validators
RE_CSV               <- "^[^,]+(,[^,]+)*$"
RE_KEY_VALUE_PAIR    <- "^[^=]+=[^=]+$"
RE_TRACK_NAME        <- "^[a-zA-Z][a-zA-Z0-9_-]*$"
RE_HEX_COLOR         <- "^#[0-9A-Fa-f]{6}$"
RE_BIGBED_TYPE       <- "^bigBed(?:\\s+([0-9]+))?(?:\\s+([+.]))?\\s*$"
RE_FILTER_COMPOSITE  <- "^[A-Za-z]+(=one)?$"
RE_SORT_ORDER        <- "^[^=]+=[+-]$"
vcf_sample_pattern   <- "^[^|,[:space:]]+(\\|[^|,[:space:]]+)?$"

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Semantic property types
###

prop_str <- function() scalar(class_character)
# relative value can be any valid string not worth checking the path is valid
# or something else. Instead handle it at runtime via append_to_uri
prop_url <- function() {
    scalar(class_character, validator = function(value) {
        has_scheme <- grepl(RE_VALID_SCHEME, value)
        if (has_scheme && !grepl(RE_HTTP_FTP_SCHEME, value))
            "must be an http/https/ftp URL or a relative path"
    })
}
prop_float <- function() scalar(class_numeric)
prop_int <- function() scalar(class_integer)
prop_literal <- function(value) scalar(class_character, choices = value)
prop_csv <- function() {
    scalar(class_character, validator = function(value) {
        if (!grepl(RE_CSV, trimws(value)))
            "must be a comma-separated list"
    })
}

on_off <- function(default = "off") {
    scalar(class_character, choices = c("on", "off"), default = default)
}

true_false <- function(default = "false") {
    scalar(class_character, choices = c("true", "false"), default = default)
}

space_sep <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) == 0L)
            "must be a non-empty whitespace-separated list"
    })
}

## "<key1>=<val1> [<key2>=<val2> ...]" — whitespace-separated key=value pairs
key_value_pairs <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) == 0L)
            return("must contain at least one <key>=<value> pair")
        if (!all(grepl(RE_KEY_VALUE_PAIR, parts)))
            "each token must be in <key>=<value> format"
    })
}

int_pair <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) != 2L)
            return("must be two whitespace-separated integers")
        vals <- suppressWarnings(as.integer(parts))
        if (any(is.na(vals)))
            "both components must be integers"
    })
}

## UCSC label-length caps are a soft spec — warn rather than hard-fail
## so parsing proceeds on hubs that violate them. S7 validators only
## signal errors (warnings raised inside validate() are swallowed) and
## S7 drops user-provided validators during class construction, so the
## cap can't live on the property. LABEL_CAPS maps property name → cap
## and check_label_caps() enforces it post-coercion in stanzas_to().
LABEL_CAPS <- list(shortLabel = 17L, longLabel = 76L)

## "<r>,<g>,<b>" — RGB triple. Named rgb_triple() to avoid shadowing
## grDevices::rgb() inside the rtracklayer namespace (used by bed.R).
rgb_triple <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), ",")[[1L]]
        if (length(parts) != 3L)
            return("must be in R,G,B format")
        vals <- suppressWarnings(as.integer(parts))
        if (any(is.na(vals)) || any(vals < 0L) || any(vals > 255L))
            "each RGB component must be an integer between 0 and 255"
    })
}

track_name <- function() {
    scalar(class_character, validator = function(value) {
        if (!grepl(RE_TRACK_NAME, value))
            "must start with a letter and contain only [a-zA-Z0-9_-]"
    })
}

track_type <- function() {
    scalar(class_character, validator = function(value) {
        base <- sub(RE_TYPE_BASE, "", trimws(value))
        if (!(base %in% TRACK_TYPES))
            paste0("unknown track type '", base, "'; expected one of: ",
                   paste(TRACK_TYPES, collapse = ", "))
    })
}

## "max:default:min"
height_range <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), ":")[[1L]]
        if (length(parts) != 3L)
            return("must be in max:default:min format")
        vals <- suppressWarnings(as.integer(parts))
        if (any(is.na(vals)) || any(vals < 0L))
            return("each component must be a non-negative integer")
        if (vals[1L] < vals[2L] || vals[2L] < vals[3L])
            "must satisfy max >= default >= min"
    })
}

## "lower:upper"
numeric_range <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), ":")[[1L]]
        if (length(parts) != 2L)
            return("must be in lower:upper format")
        vals <- suppressWarnings(as.numeric(parts))
        if (any(is.na(vals)))
            return("both components must be numeric")
        if (vals[1L] > vals[2L])
            "lower bound must not exceed upper bound"
    })
}

## "<low>" or "<low>:<high>"
score_range <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), ":")[[1L]]
        if (length(parts) < 1L || length(parts) > 2L)
            return("must be in <low> or <low>:<high> format")
        vals <- suppressWarnings(as.numeric(parts))
        if (any(is.na(vals)))
            return("components must be numeric")
        if (length(vals) == 2L && vals[1L] > vals[2L])
            "low must not exceed high"
    })
}

## bounded integer in [min, max]
bounded_int <- function(min, max) {
    scalar(class_integer, validator = function(value) {
        if (value < min || value > max)
            paste0("must be an integer between ", min, " and ", max)
    })
}

## "#RRGGBB" hex color
hex_color <- function() {
    scalar(class_character, validator = function(value) {
        if (!grepl(RE_HEX_COLOR, trimws(value)))
            "must be a hex color in #RRGGBB format"
    })
}

## "bigBed [<3-12>] [+/.]" — type-line variant for bigBed tracks
bigbed_type <- function() {
    scalar(class_character, validator = function(value) {
        m <- regmatches(value, regexec(RE_BIGBED_TYPE, trimws(value)))[[1L]]
        if (length(m) == 0L)
            return("must be 'bigBed [<3-12>] [+/.]'")
        if (nzchar(m[2L])) {
            n <- as.integer(m[2L])
            if (n < 3L || n > 12L)
                return("field count must be between 3 and 12")
        }
    })
}

## "<r,g,b> <r,g,b>" — two RGBs for colorByStrand
two_rgb <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) != 2L)
            return("must be two whitespace-separated R,G,B triples")
        for (p in parts) {
            comps <- strsplit(p, ",")[[1L]]
            if (length(comps) != 3L)
                return("each component must be in R,G,B format")
            vals <- suppressWarnings(as.integer(comps))
            if (any(is.na(vals)) || any(vals < 0L) || any(vals > 255L))
                return("each RGB component must be 0-255")
        }
    })
}

## "<path> <suffix>" — itemImagePath / itemBigImagePath
path_with_suffix <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) != 2L)
            "must be '<path> <suffix>'"
    })
}

## "<label> <url>" — downloadUrl
labeled_url <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) < 2L)
            "must be '<label> <url>'"
    })
}

## "<sampleName>" or "<sampleName>|<altName>" — VCF sample identifier
vcf_sample <- function() {
    scalar(class_character, validator = function(value) {
        if (!grepl(vcf_sample_pattern, trimws(value)))
            "must be '<sampleName>' or '<sampleName>|<altName>'"
    })
}

## comma-separated list of vcf_sample tokens
vcf_sample_list <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), ",")[[1L]]
        if (length(parts) == 0L)
            return("must contain at least one sample")
        if (!all(grepl(vcf_sample_pattern, parts)))
            "each token must be '<sampleName>' or '<sampleName>|<altName>'"
    })
}

## "<off/on/noGenome/tbNoGenome> [table1 ...]" — tableBrowser
TABLE_BROWSER_MODES <- c("off", "on", "noGenome", "tbNoGenome")
table_browser_setting <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) < 1L || !(parts[1L] %in% TABLE_BROWSER_MODES))
            paste0("first token must be one of: ",
                   paste(TABLE_BROWSER_MODES, collapse = ", "))
    })
}

## "bigChain <targetDb>" — type-line variant for bigChain tracks
bigchain_type <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) != 2L || parts[1L] != "bigChain")
            "must be 'bigChain <targetDb>'"
    })
}

## "bigWig [<min> <max>]" — type-line variant for bigWig tracks
bigwig_type <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) != 1L && length(parts) != 3L)
            return("must be 'bigWig' or 'bigWig <min> <max>'")
        if (parts[1L] != "bigWig")
            return("must start with 'bigWig'")
        if (length(parts) == 3L) {
            vals <- suppressWarnings(as.numeric(parts[2:3]))
            if (any(is.na(vals)))
                return("min and max must be numeric")
            if (vals[1L] > vals[2L])
                "min must not exceed max"
        }
    })
}

## baseColorUseSequence — accepts a small set of modes; some take extra args
BASE_COLOR_USE_SEQ_MODES <- c("extFile", "hgPcrResult", "lfExtra",
                              "nameIsSequence", "seq1Seq2", "ss", "2bit")
base_color_use_seq <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) < 1L || !(parts[1L] %in% BASE_COLOR_USE_SEQ_MODES))
            paste0("first token must be one of: ",
                   paste(BASE_COLOR_USE_SEQ_MODES, collapse = ", "))
    })
}

## baseColorDefault — five-way literal shared across multiple track types
base_color_default <- function() {
    prop_literal(c("diffBases", "diffCodons", "itemBases", "itemCodons",
              "genomicCodons"))
}

## "<gTag> <gTitle> <mTag1=mTitle1> [mTag2=mTitle2 ...]" — subGroup1..9 line.
## First two tokens are tag/title; remaining tokens are key=value members.
subgroup_def <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) < 3L)
            return("must be '<tag> <title> <mTag1=mTitle1> [mTag2=...]'")
        members <- parts[-(1:2)]
        if (!all(grepl(RE_KEY_VALUE_PAIR, members)))
            "each member token must be in <mTag>=<mTitle> format"
    })
}

## "<dim>[=one] [<dim>[=one] ...]" — filterComposite line.
## Each token is a dimension letter (A/B/C/...) optionally followed by '=one'.
filter_composite_setting <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) == 0L)
            return("must contain at least one dimension token")
        if (!all(grepl(RE_FILTER_COMPOSITE, parts)))
            "each token must be '<dim>' or '<dim>=one'"
    })
}

## "<gTag>=+/- [<gTag>=+/- ...]" — sortOrder line.
sort_order_setting <- function() {
    scalar(class_character, validator = function(value) {
        parts <- strsplit(trimws(value), RE_WHITESPACE)[[1L]]
        if (length(parts) == 0L)
            return("must contain at least one '<tag>=+/-' token")
        if (!all(grepl(RE_SORT_ORDER, parts)))
            "each token must be in '<tag>=+' or '<tag>=-' format"
    })
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Property helpers (non-scalar slot constructors and reusable groupings)
###

## class_list slot for dotted dynamic settings (filters, decorators, etc.)
dotted_list <- function() {
    new_property(class_list, default = quote(list()))
}

## "<keyword>" or integer in [min, max] — used by smoothingWindow,
## interactMultiRegion, resolution. min/max default to no bounds.
keyword_or_int <- function(keyword, min = NA_integer_, max = NA_integer_) {
    label <- if (is.na(min) && is.na(max)) "a non-negative integer"
        else if (is.na(max)) sprintf("an integer >= %d", min)
        else if (is.na(min)) sprintf("an integer <= %d", max)
        else sprintf("an integer between %d and %d", min, max)
    msg <- sprintf("must be '%s' or %s", keyword, label)
    lo <- if (is.na(min)) 0L else min
    hi <- if (is.na(max)) .Machine$integer.max else max
    nullable(scalar(class_character, validator = function(value) {
        if (identical(value, keyword)) return(NULL)
        n <- suppressWarnings(as.integer(value))
        if (is.na(n) || n < lo || n > hi) msg
    }))
}

## indelDoubleInsert / indelQueryInsert / indelPolyA — shared by bam, bigChain,
## bigPsl. Returns a named property list to splice into a class definition.
indel_props <- function() {
    list(
        indelDoubleInsert = nullable(on_off()),
        indelQueryInsert = nullable(on_off()),
        indelPolyA = nullable(on_off())
    )
}

## subGroup1..subGroup9 — composite track subgroup definitions.
## Returns a named property list to splice into the CompositeTrack definition.
subgroup_props <- function() {
    setNames(
        replicate(9L, nullable(subgroup_def()), simplify = FALSE),
        paste0("subGroup", seq_len(9L))
    )
}
