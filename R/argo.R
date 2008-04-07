## ARGO interface
## PROBLEMS:
## - Due to singletons, sessions are just view windows (w/ a tab for each view)
## - Track visibility is at the track, not the view, level

setClass("argoSession",
         representation(jref = "jobjRef", seqcache = "environment"),
         contains = "browserSession")

setMethod("initialize", "argoSession",
          function(.Object, name = "R Session",
                   jar = system.file("java", "argo.jar",
                     package = "rtracklayer"))
          {
            .jinit(jar)
            frame <- .jcall("calhoun/gebo/ui/DesktopFrame",
                            "Lcalhoun/gebo/ui/DesktopFrame;", "getInstance")
            if (is.null(frame)) {
### HACK: block garbage on System.out
### Really need to subclass OutputStream with no-op
              out <- .jnew("java/io/PrintStream", tempfile())
              .jcall("java/lang/System", "V", "setOut", out)
              
              .jcall("calhoun/gebo/ui/DesktopFrame", "V", "main",
                     .jarray(character(0)), FALSE)
              frame <- .jcall("calhoun/gebo/ui/DesktopFrame",
                              "Lcalhoun/gebo/ui/DesktopFrame;", "getInstance")
              window <- .jcall(frame, "Ljavax/swing/JInternalFrame;",
                               "getSelectedFrame")
            } else {
              window <- .jcall("calhoun/gebo/action/SequenceViewAction",
                               "Lcalhoun/gebo/ui/SegmentViewWindow;", "open")
            }
            .jcall(window, "V", "setTitle", name)
            .Object@jref <- window
            .Object
          })

setMethod("browserViews", "argoSession",
          function(object)
          {
            maps <- .jcall(object@jref,
                           "[Lcalhoun/gebo/ui/featuremap/FeatureMap;",
                           "getFeatureMaps")
            lapply(maps,
                   function(map) new("argoView", jref = map, session = object))
          })

setMethod("activeView", "argoSession",
          function(object)
          {
            map <- .jcall(object@jref,
                          "Lcalhoun/gebo/ui/featuremap/FeatureMap;",
                          "getActiveFeatureMap")
            if (!is.null(map))
              new("argoView", jref = map, session = object)
            else NULL
          })

setReplaceMethod("activeView", "argoSession",
                 function(object, value)
                 {
                   active <- activeView(object)
                   tabbed <- .jcall(active@jref, "Ljavax/swing/JTabbedPane;",
                                    "getTabbedPane")
                   .jcall(tabbed, "V", "setSelectedComponent",
                          .jcast(value@jref, "java/awt/Component"))
                   object
                 })

setMethod("tracks", "argoSession",
          function(object, segment = NULL, visible = FALSE)
          {
            tracks <- argoTracks(object, segment, visible)
            sapply(tracks, .jcall, "S", "getId")
          })

setMethod("trackSet", "argoSession",
          function(object, segment = genomeSegment(object), name)
          {
            registry <- .jcall("calhoun/gebo/db/TrackManagerRegistry",
                               "Lcalhoun/gebo/db/TrackManagerRegistry;",
                               "getGlobal")
            managers <- .jcall(registry, "[Lcalhoun/gebo/api/TrackManager;",
                               "getTrackManagers")
            for (manager in managers) {
              print(manager)
              track <- .jcall(manager, "Lcalhoun/gebo/model/FeatureTrack;",
                              "getTrack", name)
              if (!is.null(track))
                break
            }
            jsegment <- .jcast(argoJobj(segment), "calhoun/gebo/model/Segment")
            features <- .jcall(track, "[Lcalhoun/gebo/model/Feature;",
                               "getFeatures", jsegment)
            getFeatureRow <- function(feature)
              {
                featStart <- .jcall(feature, "I", "getStart")
                featEnd <- .jcall(feature, "I", "getStop")
                strand <- .jcall(feature, "Lcalhoun/gebo/model/Strand;",
                                 "getStrand")
                featStrand <- .jcall(strand, "S", "toString")
                data.frame(featStart = featStart, featEnd = featEnd,
                           featStrand = featStrand)
              }
            featureData <- do.call("rbind", lapply(features, getFeatureRow))
            featureData <- cbind(featureData, featChrom = segment@chrom)
            new("trackSet", featureData = featureData,
                dataVals = rep(NA, nrow(featureData)),
                genome = segment@genome)
          })

setMethod("layTrack", c("argoSession", "trackSet"),
          function(object, track, name, view)
          {
            argoTrack(track, object, name, view)
            if (view)
              browserView(object, genomeSegment(track))
            object
          })

setMethod("laySequence", c("argoSession", "character"),
          function(object, sequence, name, label = name)
          {
            jsequence <- .jnew("calhoun/gebo/model/Sequence", sequence, name,
                               label)
            registry <- .jcall("calhoun/gebo/db/SequenceManagerRegistry",
                               "Lcalhoun/gebo/db/SequenceManagerRegistry;",
                               "getGlobal")
            manager <- .jcall(jsequence, "Lcalhoun/gebo/api/SequenceManager;",
                              "getSequenceManager")
            .jcall(registry, "V", "add", manager)
            object
          })

setMethod("genomeSequence", "argoSession",
          function(object, segment)
          {
            jsegment <- argoJobj(segment, object)
            jsequence <- .jcall(jsegment, "Lcalhoun/gebo/model/Sequence;",
                                "getSequence")
            jmanager <- .jcall(jsequence, "Lcalhoun/gebo/api/SequenceManager;",
                               "getSequenceManager")
            BString(.jcall(jmanager, "S", "getRaw", jsequence))
          })

setMethod("close", "argoSession",
          function(con) {
            frame <- .jcall("calhoun/gebo/ui/DesktopFrame",
                            "Lcalhoun/gebo/ui/DesktopFrame;", "getInstance")
            desktop <- .jcall(frame, "Ljavax/swing/JDesktopPane;",
                              "getDesktopPane")
            .jcall(desktop, "V", "remove", con@jref)
          })

setClass("argoView", representation(jref = "jobjRef"),
         contains = "browserView")

setMethod("browserView", "argoSession",
          function(object, segment, track = tracks(object, segment, TRUE), ...)
          {
            if (length(list(...))) {
              if (is.null(segment))
                segment <- genomeSegment(...)
              else segment <- genomeSegment(..., segment = segment)
            }
            jsegment <- argoJobj(segment, object)
            .jcall(object@jref, "V", "addFeatureMap",
                   .jcast(jsegment, "calhoun/gebo/model/Segment"),
                   .jarray(list(), "calhoun/gebo/model/Feature"))
            view <- activeView(object)
            tracks(view) <- track
            view
          })

setMethod("activeView", "argoView",
          function(object)
          {
            tabbed <- .jcall(active, "Ljavax/swing/JTabbedPane;",
                             "getTabbedPane")
            sel <- .jcall(tabbed, "Ljava/awt/Component;",
                          "getSelectedComponent")
            .jequals(sel, object@jref)
          })

setMethod("close", "argoView",
          function(con)
          {
            .jcall(con@session@jref, "V", "removeFeatureMap",
                   .jcast(con@jref, "calhoun/gebo/ui/featuremap/FeatureMap"))
          })

setMethod("tracks", "argoView",
          function(object)
          {
            tracks <- .jcall(object@jref, "Ljava/util/List;",
                             "getVisibleTracks")
            tracks <- .jcall(tracks, "[Ljava/lang/Object;", "toArray",
                             .jarray(list(), "java/lang/Object"))
            sapply(tracks, .jcall, "S", "getId")
          })
setReplaceMethod("tracks", "argoView",
                 function(object, value)
                 {
                   tracks <- argoTracks(object)
                   ids <- lapply(tracks, .jcall, "S", "getId")
                   visible <- ids %in% value
                   for (i in seq_along(ids))
                     .jcall(tracks[[i]], "V", "setVisible", visible[i])
                   msg <- .jfield("calhoun/gebo/ui/featuremap/FeatureMap",
                                  "S", "ALL_CHANGE", convert = FALSE)
                   .jcall("calhoun/gebo/ui/featuremap/FeatureMap", "V",
                          "notify", argoSequence(object),
                          .jcast(msg, "java/lang/Object"))
                   object
                 })

setMethod("genomeSegment", "argoView",
          function(object)
          {
            segment <- .jcall(object@jref, "Lcalhoun/gebo/model/Segment;",
                              "getSegment")
            start <- .jcall(segment, "I", "getStart")
            end <- .jcall(segment, "I", "getStop")
            sequence <- .jcall(segment, "Lcalhoun/gebo/model/Sequence;",
                               "getSequence")
            genome <- character(0)
            id <- .jcall(sequence, "S", "getId")
            group <- .jcall(sequence, "Lcalhoun/gebo/model/SequenceGroup;",
                            "getPrimarySequenceGroup")
            if (!is.null(group)) {
              genome <- .jcall(group, "S", "getId")
              chrom <- id
            } else {
              idSplit <- strsplit(id, " ")[[1]]
              chrom <- tail(idSplit, 1)
              if (length(idSplit) > 1)
                genome <- idSplit[1]
            }
            genomeSegment(genome = genome, chrom = chrom,
                          start = start, end = end)
          })

setReplaceMethod("genomeSegment", "argoView",
                 function(object, value)
                 {
                   # we can't change the sequence, only start and stop
                   segment <- .jcall(object@jref, "Lcalhoun/gebo/model/Segment",
                                     "getSegment")
                   .jcall(segment, "V", "setStart", as.integer(value@start))
                   .jcall(segment, "V", "setStop", as.integer(value@end))
                   object
                 })

# something like this
setMethod("selectedFeatures", "argoView",
          function(object)
          {
            tracks <- argoTracks(object)
            selected <- function(track)
              {
                features <- .jcall(track, "[Lcalhoun/gebo/model/Feature;",
                                   "getFeatures", sequence)
                lapply(features, .jcall, "Z", "isSelected")
              }
            sel <- lapply(tracks, selected)
            names(sel) <- sapply(tracks, .jcall, "S", "getId")
            sel
          })

# utilities

setGeneric("argoJobj", function(object, ...) standardGeneric("argoJobj"))

## resolve to a calhoun.gebo.model.Segment reference
setMethod("argoJobj", "genomeSegment",
          function(object, session)
          {
            ## ID is 'genome chrom:start-end'
            #browser()
            id <- object@chrom
            if (length(object@genome))
              id <- paste(object@genome, id)

            jseq <- NULL
            ## lookup in our custom sequence manager (or an R env)
            if (!missing(session))
              jseq <- session@seqcache[[id]]
            ##manager <- .jcall("org/bioc/rtracklayer/RSequenceManager",
            ##                  "Lorg/bioc/rtracklayer/RSequenceManager;",
            ##                  "getInstance")
            ## jseq <- .jcall(manager, "Lcalhoun/gebo/model/Sequence;",
            ##                "getSequence", id)
            if (is.null(jseq)) {
              ## if a sequence does not exist, lookup in UCSC
              ucsc <- .jcall("calhoun/gebo/db/UCSCDataManager",
                             "Lcalhoun/gebo/db/UCSCDataManager;",
                             "getInstance")
              jseq <- .jcall(ucsc, "Lcalhoun/gebo/model/Sequence;",
                             "getSequence", id)
              if (is.null(jseq))
                stop("Could not find sequence with ID: ", id)
### TODO: and add to our manager
              ##.jcall(manager, "V", "addSequence", sequence, id)
            }
            if (!missing(session))
              session@seqcache[[id]] <- jseq
            
            .jnew("calhoun/gebo/model/SimpleSegment", jseq,
                  as.integer(object@start), as.integer(object@end))
          })

setGeneric("argoSequence",
           function(object, ...) standardGeneric("argoSequence"))

setMethod("argoSequence", "argoView",
          function(object)
          {
            segment <- .jcall(object@jref, "Lcalhoun/gebo/model/Segment;",
                              "getSegment")
            .jcall(segment, "Lcalhoun/gebo/model/Sequence;", "getSequence")
          })

setGeneric("argoTracks", function(object, ...) standardGeneric("argoTracks"))

setMethod("argoTracks", "argoView",
          function(object)
          {
            sequence <- argoSequence(object)
            .jcall(sequence, "[Lcalhoun/gebo/model/FeatureTrack;", "getTracks")
          })

setMethod("argoTracks", "argoSession",
          function(object, segment = NULL, visible = FALSE)
          {
            registry <- .jcall("calhoun/gebo/db/TrackManagerRegistry",
                               "Lcalhoun/gebo/db/TrackManagerRegistry;",
                               "getGlobal")
            ## method appears to be broken
            ##tracks <- .jcall(registry, "[Lcalhoun/gebo/model/FeatureTrack;",
            ##                 "getTracks")
            managers <- .jcall(registry, "[Lcalhoun/gebo/api/TrackManager;",
                               "getTrackManagers")
            if (is.null(segment))
              trackList <- lapply(managers, .jcall,
                                  "[Lcalhoun/gebo/model/FeatureTrack;",
                                  "getTracks")
            else {
              sequence <- .jcall(argoJobj(segment, object),
                                 "Lcalhoun/gebo/model/Sequence;", "getSequence")
              trackList <- lapply(managers, .jcall,
                                  "[Lcalhoun/gebo/model/FeatureTrack;",
                                  "getTracks", sequence)
            }
            tracks <- unlist(trackList)
            if (visible) {
              visibleTracks <- sapply(tracks, .jcall, "Z", "isVisible")
              tracks <- tracks[visibleTracks]
            }
            tracks
          })

setGeneric("argoTrack", function(object, ...) standardGeneric("argoTrack"))

setMethod("argoTrack", "trackSet",
          function(object, session, name, visible = FALSE)
          {
            manager <- .jcall("calhoun/gebo/db/DynamicTrackManager",
                              "Lcalhoun/gebo/db/DynamicTrackManager;",
                              "getInstance")
            df <- trackData(object)
            segment <- genomeSegment(object)
            sequence <- .jcall(argoJobj(segment, session),
                               "Lcalhoun/gebo/model/Sequence;", "getSequence")
            if (is.null(df$color)) {
              black <- .jfield("java/awt/Color", "Ljava/awt/Color;", "BLACK")
              color <- rep(list(black), nrow(df))
            } else {
              jcolor <- function(rgb)
                .jnew("java/awt/Color", rgb[1], rgb[2], rgb[3])
              color <- apply(col2rgb(df$color), 2, jcolor)
            }
            if (is.null(df$featStrand))
              df$featStrand <- "+"
            for (i in seq_len(nrow(df))) {
### FIXME: how to add data values?
              strand <- .jcall("calhoun/gebo/model/Strand",
                               "Lcalhoun/gebo/model/Strand;", "parseStrand",
                               df$featStrand[i]) # FIXME: inefficient
              segment <- .jnew("calhoun/gebo/model/SimpleSegment", sequence,
                               strand, df$featStart[i], df$featEnd[i])
              segment <- .jcast(segment, "calhoun/gebo/model/Segment")
              feature <- .jcall(manager, "Lcalhoun/gebo/model/Feature;",
                                "createFeature", segment, color[[i]], name)
            }
                                        # make track visible
            track <- .jcall(feature, "Lcalhoun/gebo/model/Track;",
                            "getTrack")
            .jcall(track, "V", "setVisible", visible)
            .jcall("calhoun/gebo/ui/featuremap/FeatureMap", "V", "notify",
                   sequence, .jcast(manager, "java/lang/Object"))
            track
          })
