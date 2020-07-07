### =========================================================================
### IGV Data Server Access
### -------------------------------------------------------------------------

readIGVRegistry <- function(url) {
    do.call(rbind, lapply(readLines(url), readIGVDataFile))
}

handleCategory <- function(x) {
    paths <- unlist(getNodeSet(x, "Resource/@path"))
    names(paths) <- unlist(getNodeSet(x, "Resource/@name"))
    name <- factor(unname(xmlAttrs(x)["name"]))
    categories <- lapply(getNodeSet(x, "Category"), handleCategory)
    children <- do.call(rbind, categories)
    children$parenet <- name
    resources <- BiocFileList(lapply(paths, FileForFormat))
    df <- DataFrame(category=name, parent=NA, name=names(paths), resources)
    rbind(df, children)
}

readIGVDataFile <- function(url) {
    lines <- readLines(url)
    doc <- xmlTreeParse(lines, asText=TRUE, useInternalNodes=TRUE)
    handleCategory(getNodeSet(doc, "/Global")[[1L]])
}
