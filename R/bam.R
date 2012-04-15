### =========================================================================
### BAM support (wrappers around Rsamtools)
### -------------------------------------------------------------------------

setMethod("import", "BamFile",
          function(con, format, text, use.names = FALSE,
                   param = ScanBamParam(...), ...)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            readBamGappedAlignments(con, use.names = use.names, param = param)
          })
