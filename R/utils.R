pasteCollapse <- function(x, collapse = ",", sentinel = "\r") {
  s <- start(PartitioningByWidth(x))
  x_flat <- unlist(x, use.names = FALSE)
  x_flat[s] <- paste(sentinel, x_flat[s], sep = "")
  big_list <- paste(x_flat, collapse = collapse)
  gsub_str <- paste(collapse, sentinel, sep="")
  big_list <- gsub(gsub_str, sentinel, big_list, fixed = TRUE)
  big_list <- substring(big_list, nchar(sentinel))
  unlist(strsplit(big_list, sentinel, fixed = TRUE))
}
