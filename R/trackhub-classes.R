### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Cross-object reference validation
###

Refer := new_class(
    properties = list(
        valid = new_property(class_character, default = character(0L)),
        child = scalar(class_any),
        child_key = prop_str(),
        context = nullable(prop_str())
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Hub & Genome classes
###

Hub := new_class(
    properties = list(
        hub = prop_str(),
        shortLabel = prop_str(),
        longLabel = prop_str(),
        genomesFile = scalar(class_character, default = GENOME_FILE_NAME),
        email = prop_str(),
        descriptionUrl = nullable(prop_url()),
        useOneFile = nullable(on_off())
    )
)

Genome := new_class(
    properties = list(
        genome = prop_str(),
        trackDb = scalar(class_character, default = TRACKDB_FILE_NAME),
        metaDb = nullable(scalar(class_character, default = METADB_FILE_NAME)),
        metaTab = nullable(scalar(class_character, default = METATAB_FILE_NAME))
    )
)

Assembly := new_class(
    Genome,
    properties = list(
        twoBitPath = prop_url(),
        groups = nullable(scalar(class_character, default = GROUP_FILE_NAME)),
        description = nullable(prop_str()),
        organism = nullable(prop_str()),
        defaultPos = nullable(prop_str()),
        orderKey = nullable(prop_str()),
        htmlPath = nullable(prop_url())
    )
)

Genome_OR_Assembly <- Genome | Assembly

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Track classes
###

CommonSettings := new_class(
    properties = list(
        track = track_name(),
        ## caps enforced via LABEL_CAPS in check_label_caps()
        shortLabel = prop_str(),
        longLabel = prop_str(),
        visibility = nullable(
            scalar(class_character,
                   choices = c("hide", "dense", "squish", "pack", "full"),
                   default = "hide")
        ),
        html = nullable(prop_url())
    )
)

CommonOptionalSettings := new_class(
    parent = CommonSettings,
    properties = list(
        color = nullable(rgb_triple()),
        priority = nullable(prop_float()),
        altColor = nullable(rgb_triple()),
        boxedCfg = nullable(on_off()),
        chromosomes = nullable(prop_csv()),
        darkerLabels = nullable(prop_literal("on")),
        dataVersion = nullable(prop_str()),
        directUrl = nullable(prop_url()),
        downloadUrl = nullable(labeled_url()),
        iframeUrl = nullable(prop_url()),
        iframeOptions = nullable(prop_str()),
        mouseOver = nullable(prop_str()),
        mouseOverField = nullable(prop_str()),
        multiRegionsBedUrl = nullable(prop_url()),
        otherDb = nullable(prop_str()),
        otherTwoBitUrl = nullable(prop_url()),
        pennantIcon = nullable(prop_str()),
        tableBrowser = nullable(table_browser_setting()),
        url = nullable(prop_url()),
        urlLabel = nullable(prop_str()),
        urls = nullable(prop_str()),
        skipEmptyFields = nullable(prop_literal("on")),
        skipFields = nullable(prop_csv()),
        sepFields = nullable(prop_csv())
    )
)

## Track: base class for any stanza that isn't a container (superTrack,
## compositeTrack, multiWig overlay container, or view). Cross-object
## reference validation for `group` and `parent` lives on TrackHub — a
## Track in isolation doesn't know about its hub's groups or sibling
## tracks, so it can't check those references itself.
Track := new_class(
    CommonOptionalSettings,
    properties = list(
        type = nullable(track_type()),
        bigDataUrl = nullable(prop_url()),
        meta = nullable(prop_str()),
        parent = nullable(prop_str()),
        group = nullable(prop_str()),
        ## subgroup membership for subtracks of a composite
        subGroups = nullable(key_value_pairs()),
        ## "Item or region" settings broadly available across data tracks
        maxWindowCoverage = nullable(prop_int()),
        maxWindowToDraw = nullable(prop_int())
    )
)

SuperTrack := new_class(
    CommonOptionalSettings,
    properties = list(
        superTrack = prop_literal(c("on", "on show"))
    )
)

CompositeTrack := new_class(
    CommonOptionalSettings,
    properties = c(
        list(
            ## "on" for a regular composite, "faceted" for a faceted composite
            compositeTrack = prop_literal(c("on", "faceted")),
            parent = nullable(prop_str()),
            allButtonPair = nullable(on_off()),
            centerLabelsDense = nullable(on_off()),
            dragAndDrop = nullable(prop_literal("subTracks")),
            ## hide-empty-subtracks family
            hideEmptySubtracks = nullable(on_off()),
            hideEmptySubtracksMultiBedUrl = nullable(prop_url()),
            hideEmptySubtracksSourcesUrl = nullable(prop_url()),
            hideEmptySubtracksLabel = nullable(prop_str())
        ),
        ## subGroup1..subGroup9
        subgroup_props(),
        list(
            ## subgroup display/selection settings
            dimensions = nullable(key_value_pairs()),
            filterComposite = nullable(filter_composite_setting()),
            dimensionAchecked = nullable(space_sep()),
            dimensionBchecked = nullable(space_sep()),
            dimensionCchecked = nullable(space_sep()),
            sortOrder = nullable(sort_order_setting()),
            ## faceted composite settings (compositeTrack = "faceted")
            metaDataUrl = nullable(prop_url()),
            primaryKey = nullable(prop_str()),
            maxCheckBoxes = nullable(prop_int()),
            dataTypes = nullable(space_sep())
        )
    )
)

OverlayTrack := new_class(
    CommonOptionalSettings,
    properties = list(
        container = prop_literal("multiWig"),
        parent = nullable(prop_str()),
        aggregate = prop_literal(c("transparentOverlay", "stacked",
                              "solidOverlay", "none")),
        showSubtrackColorOnUi = nullable(on_off())
    )
)

## A "view" stanza inside a multi-view composite. Discriminated by the
## presence of a `view` key (the view tag). Children reference this stanza
## via `parent <viewName>`.
View := new_class(
    CommonOptionalSettings,
    properties = list(
        view = prop_str(),
        parent = nullable(prop_str()),
        viewUi = nullable(on_off()),
        configurable = nullable(on_off())
    )
)

BamTrack := new_class(
    Track,
    properties = c(
        list(
            refUrl = nullable(prop_url()),
            bigDataIndex = nullable(prop_url()),
            bamColorMode = nullable(prop_literal(c("strand", "gray", "tag", "off"))),
            bamGrayMode = nullable(prop_literal(c("aliQual", "baseQual",
                                             "unpaired"))),
            aliQualRange = nullable(numeric_range()),
            baseQualRange = nullable(numeric_range()),
            bamColorTag = nullable(prop_str()),
            noColorTag = nullable(prop_literal(".")),
            bamSkipPrintQualScore = nullable(prop_literal("."))
        ),
        indel_props(),
        list(
            minAliQual = nullable(prop_int()),
            pairEndsByName = nullable(prop_literal(".")),
            pairSearchRange = nullable(prop_int()),
            showNames = nullable(on_off()),
            doWiggle = nullable(prop_literal("on"))
        )
    )
)

BigBarChartTrack := new_class(
    Track,
    properties = list(
        barChartBars = nullable(space_sep()),
        barChartColors = nullable(space_sep()),
        barChartLabel = nullable(prop_str()),
        barChartMaxSize = nullable(prop_literal(c("small", "medium", "large"))),
        barChartSizeWindows = nullable(int_pair()),
        barChartStretchToItem = nullable(on_off()),
        barChartFacets = nullable(on_off()),
        barChartStatsUrl = nullable(prop_url()),
        singleCellColumnNames = nullable(on_off()),
        barChartMerge = nullable(on_off()),
        barChartMetric = nullable(prop_str()),
        barChartUnit = nullable(prop_str()),
        barChartCategoryUrl = nullable(prop_url()),
        barChartSampleUrl = nullable(prop_url()),
        barChartBarMinPadding = nullable(prop_float()),
        barChartBarMinWidth = nullable(prop_float()),
        maxLimit = nullable(prop_float()),
        labelFields = nullable(prop_csv()),
        defaultLabelFields = nullable(prop_csv())
    )
)

BigBedTrack := new_class(
    Track,
    properties = list(
        ## type override: bigBed accepts an extended type spec
        type = nullable(bigbed_type()),
        ## color/display
        itemRgb = nullable(on_off()),
        colorByStrand = nullable(two_rgb()),
        denseCoverage = nullable(prop_int()),
        labelOnFeature = nullable(on_off()),
        exonArrows = nullable(on_off()),
        exonNumbers = nullable(on_off()),
        minGrayLevel = nullable(bounded_int(1L, 9L)),
        spectrum = nullable(prop_literal("on")),
        style = nullable(prop_literal("heatmap")),
        thickDrawItem = nullable(on_off()),
        ## scoring
        scoreFilter = nullable(score_range()),
        noScoreFilter = nullable(prop_literal("on")),
        ## detail tables
        extraDetailsTable = nullable(prop_url()),
        extraTableFields = nullable(prop_csv()),
        detailsStaticTable = nullable(prop_url()),
        detailsDynamicTable = nullable(prop_csv()),
        ## item limits
        maxItems = nullable(prop_int()),
        ## search & labels
        searchIndex = nullable(prop_str()),
        searchTrix = nullable(prop_url()),
        labelFields = nullable(prop_csv()),
        defaultLabelFields = nullable(prop_csv()),
        labelSeparator = nullable(prop_str()),
        ## highlight
        highlightColor = nullable(hex_color()),
        ## dynamic dotted settings — keyed by field name
        ## filter.<f>, filterText.<f>, filterValues.<f>, filterLabel.<f>
        filters = dotted_list(),
        ## highlight.<f>, highlightText.<f>, highlightValues.<f>
        highlights = dotted_list(),
        ## decorator.<...>
        decorators = dotted_list(),
        ## less frequent
        bedNameLabel = nullable(prop_str()),
        exonArrowsDense = nullable(on_off()),
        itemImagePath = nullable(path_with_suffix()),
        mergeSpannedItems = nullable(on_off()),
        linkIdInName = nullable(prop_literal("on")),
        nextExonText = nullable(prop_str()),
        scoreLabel = nullable(prop_str()),
        showTopScorers = nullable(prop_int())
    )
)

BigChainTrack := new_class(
    Track,
    properties = c(
        list(
            ## type override: bigChain takes a target database name
            type = nullable(bigchain_type()),
            linkDataUrl = nullable(prop_url()),
            baseColorUseSequence = nullable(base_color_use_seq()),
            baseColorDefault = nullable(base_color_default())
        ),
        indel_props()
    )
)

BigGenePredTrack := new_class(
    Track,
    properties = list(
        baseColorDefault = nullable(base_color_default()),
        labelFields = nullable(prop_csv()),
        defaultLabelFields = nullable(prop_csv()),
        labelSeparator = nullable(prop_str()),
        ## decorator.<...> — dotted dynamic, keyed by suffix
        decorators = dotted_list()
    )
)

BigInteractTrack := new_class(
    Track,
    properties = list(
        interactDirectional = nullable(prop_literal(c("true", "offsetSource",
                                                 "offsetTarget",
                                                 "clusterSource",
                                                 "clusterTarget"))),
        interactUp = nullable(true_false()),
        detailsBoxesEnabled = nullable(true_false()),
        interactMultiRegion = keyword_or_int("true"),
        endsVisible = nullable(prop_literal("two")),
        maxHeightPixels = nullable(height_range()),
        scoreMin = nullable(prop_int()),
        spectrum = nullable(prop_literal("on"))
    )
)

BigLollyTrack := new_class(
    Track,
    properties = list(
        noStems = nullable(on_off()),
        lollySizeField = nullable(prop_int()),
        lollyMaxSize = nullable(prop_int()),
        lollyField = nullable(prop_int()),
        ## yAxisLabel.<integer> — dotted dynamic, keyed by suffix
        yAxisLabels = dotted_list(),
        yAxisNumLabels = nullable(on_off())
    )
)

BigMafTrack := new_class(
    Track,
    properties = list(
        speciesOrder = nullable(space_sep()),
        speciesLabels = nullable(key_value_pairs()),
        frames = nullable(prop_url()),
        summary = nullable(prop_str())
    )
)

BigNarrowPeakTrack := new_class(
    Track,
    properties = list(
        pValueFilter = nullable(score_range()),
        qValueFilter = nullable(score_range()),
        signalFilter = nullable(score_range())
    )
)

BigPslTrack := new_class(
    Track,
    properties = c(
        list(
            baseColorUseCds = nullable(prop_literal("given")),
            baseColorUseSequence = nullable(base_color_use_seq()),
            baseColorDefault = nullable(base_color_default()),
            showDiffBasesAllScales = nullable(prop_literal("on"))
        ),
        indel_props(),
        list(
            pslSequence = nullable(prop_literal(c("no", "all", "different"))),
            showCdsAllScales = nullable(prop_literal("on")),
            showCdsMaxZoom = nullable(prop_float()),
            showDiffBasesMaxZoom = nullable(prop_float()),
            labelFields = nullable(prop_csv()),
            defaultLabelFields = nullable(prop_csv()),
            labelSeparator = nullable(prop_str()),
            ## decorator.<...> — dotted dynamic, keyed by suffix
            decorators = dotted_list()
        )
    )
)

BigWigTrack := new_class(
    Track,
    properties = list(
        type = nullable(bigwig_type()),
        autoScale = nullable(prop_literal(c("off", "on", "group"))),
        maxHeightPixels = nullable(height_range()),
        viewLimits = nullable(numeric_range()),
        viewLimitsMax = nullable(numeric_range()),
        alwaysZero = nullable(on_off()),
        graphTypeDefault = nullable(prop_literal("points")),
        maxWindowToQuery = nullable(prop_int()),
        negateValues = nullable(prop_literal("on")),
        setColorWith = nullable(prop_url()),
        smoothingWindow = keyword_or_int("off", 1L, 16L),
        transformFunc = nullable(prop_literal(c("NONE", "LOG"))),
        logo = nullable(prop_literal("on")),
        logoMaf = nullable(prop_url()),
        windowingFunction = nullable(prop_literal(c("mean", "mean+whiskers",
                                               "maximum", "minimum"))),
        yLineMark = nullable(prop_float()),
        yLineOnOff = nullable(on_off()),
        gridDefault = nullable(prop_literal("on"))
    )
)

HalSnakeTrack := new_class(
    Track,
    properties = list(
        showSnpWidth = nullable(prop_int()),
        otherSpecies = nullable(prop_str())
    )
)

HicTrack := new_class(
    Track,
    properties = list(
        drawMode = nullable(prop_literal(c("triangle", "square", "arc"))),
        normalization = nullable(prop_literal(c("NONE", "VC", "VC_SQRT", "KR"))),
        resolution = keyword_or_int("Auto"),
        saturationScore = nullable(prop_float()),
        autoScale = nullable(prop_literal(c("off", "on", "group"))),
        hicDistanceMin = nullable(prop_int()),
        hicDistanceMax = nullable(prop_int()),
        hicArcLimit = nullable(prop_int()),
        hicArcLimitEnabled = nullable(true_false())
    )
)

VcfTabixTrack := new_class(
    Track,
    properties = list(
        bigDataIndex = nullable(prop_url()),
        ## haplotype clustering
        hapClusterEnabled = nullable(true_false()),
        ## "centerWeighted" | "fileOrder" | "treeFile <url>"
        hapClusterMethod = nullable(
            scalar(class_character, validator = function(value) {
                parts <- strsplit(trimws(value), "\\s+")[[1L]]
                mode <- parts[1L]
                if (mode %in% c("centerWeighted", "fileOrder")) {
                    if (length(parts) != 1L)
                        return(paste0("'", mode, "' takes no arguments"))
                    return(NULL)
                }
                if (identical(mode, "treeFile")) {
                    if (length(parts) != 2L)
                        return("'treeFile' must be followed by a single url")
                    return(NULL)
                }
                "must be 'centerWeighted', 'fileOrder', or 'treeFile <url>'"
            })
        ),
        hapClusterColorBy = nullable(prop_literal(c("altOnly", "function",
                                               "refAlt", "base"))),
        geneTrack = nullable(prop_str()),
        hapClusterTreeAngle = nullable(prop_literal(c("triangle", "rectangle"))),
        hapClusterHeight = nullable(prop_int()),
        ## quality / frequency
        applyMinQual = nullable(true_false()),
        minQual = nullable(prop_float()),
        minFreq = nullable(prop_float()),
        vcfDoFilter = nullable(on_off()),
        vcfDoQual = nullable(on_off()),
        vcfDoMaf = nullable(on_off())
    )
)

VcfPhasedTrioTrack := new_class(
    Track,
    properties = list(
        vcfChildSample = nullable(vcf_sample()),
        bigDataIndex = nullable(prop_url()),
        vcfParentSamples = nullable(vcf_sample_list()),
        vcfUseAltSampleNames = nullable(on_off()),
        geneTrack = nullable(prop_str()),
        vcfDoFilter = nullable(on_off()),
        vcfDoQual = nullable(on_off()),
        vcfDoMaf = nullable(on_off())
    )
)

AnyTrack <- Track | SuperTrack | CompositeTrack | OverlayTrack | View |
    BamTrack | BigBarChartTrack | BigBedTrack | BigChainTrack |
    BigGenePredTrack | BigInteractTrack | BigLollyTrack | BigMafTrack |
    BigNarrowPeakTrack | BigPslTrack | BigWigTrack |
    HalSnakeTrack | HicTrack | VcfTabixTrack | VcfPhasedTrioTrack

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Group & MetaData classes
###

Group := new_class(
    properties = list(
        name = prop_str(),
        label = nullable(prop_str()),
        priority = nullable(prop_float()),
        defaultIsClosed = nullable(scalar(class_integer, choices = c(0L, 1L)))
    )
)

MetaData := new_class(
    properties = list(
        table = new_property(class_data.frame,
                             default = quote(data.frame()))
    ),
    validator = function(self) {
        if (nrow(self@table) > 0L && !("meta" %in% colnames(self@table)))
            "table must contain a 'meta' column"
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TrackHub container
###
### Instead of making TrackHub hierarchical, we keep the object flat
### and give the illusion of hub$genome$tracks via accessors.

TrackHub := new_class(
    properties = list(
        hub = scalar(Hub),
        genomes = list_of(Genome_OR_Assembly, named = TRUE),
        tracks = list_of(AnyTrack, named = TRUE),
        metadata = nullable(scalar(MetaData)),
        groups = nullable(list_of(Group, named = TRUE))
    ),
    validator = function(self) {
        ## Only data-track stanzas carry a `group` slot — containers
        ## (SuperTrack, CompositeTrack, OverlayTrack, View) do not.
        group_names <- if (length(self@groups) > 0L)
            vapply(self@groups, function(g) g@name, character(1L))
        else character(0L)
        msgs <- unlist(lapply(self@tracks, function(track) {
            if (!S7_inherits(track, Track)) return(NULL)
            is_valid(Refer(valid = group_names,
                           child = track,
                           child_key = "group",
                           context = paste0("track '", track@track, "'")))
        }))
        if (length(msgs) > 0L) msgs
    }
)
