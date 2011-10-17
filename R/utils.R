normURI <- function(x) {
  if (!isSingleString(uri))
    stop("URI must be a single, non-NA string")
  uri <- parseURI(x)
  if (!nzchar(uri$scheme))
    x <- paste("file://", file_path_as_absolute(x), sep = "")
  x
}

pasteCollapse <- function(x, collapse = ",", sentinel = "\r") {
  s <- start(PartitioningByWidth(x))
  x_flat <- unlist(x, use.names = FALSE)
  x_flat[s] <- paste(sentinel, x_flat[s], sep = "")
  big_list <- paste(v, collapse = collapse)
  gsub_str <- paste(concat_str, sentinel, sep="")
  big_list <- gsub(gsub_str, sentinel, big_list, fixed = TRUE)
  big_list <- substring(big_list, nchar(sentinel))
  unlist(strsplit(big_list, sentinel, fixed = TRUE))
}
