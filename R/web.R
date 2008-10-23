htmlErrorHandler <- function(msg, code, domain, line, col, level, filename) {
  if (level > 2)
    stop("Failed to Parse HTML [", line, ":", col, "]: ", msg)
}

htmlParse <- function(str)
  suppressWarnings(htmlTreeParse(str, asText = TRUE, useInternalNodes = TRUE,
                                 error = htmlErrorHandler))

httpGet <- function(url, .form = list(), .parse = TRUE, ...) {
  if (length(.form) == 0)
    out <- getURL(url, ...)
  else out <- getForm(url, .params = .form, .opts = list(...))
  if (.parse)
    htmlParse(out)
  else out
}

httpPost <- function(url, .form = list(), .parse = TRUE, ...) {
  form <- postForm(url, .params = .form, .opts = list(...))
  if (.parse)
    htmlParse(form)
  else form
}

httpShow <- function(url, .form = list(), ...)
  browseURL(urlForm(url, .form, ...))

urlForm <- function(url, .form = list(), ...) {
  values <- urlEncode(as.character(c(.form, ...)))
  query <- paste(names(.form), values, sep = "=", collapse = "&")
  paste(url, query, sep = "?")
}

# differs from URLencode (vectorized, only encodes specified chars)
urlEncode <- function(str, chars = "-a-zA-Z0-9$_.+!*'(),", keep = TRUE)
{
  bad <- chars
  if (keep)
    bad <- gsub(paste("[", chars, "]", sep=""), "", str)
  bad <- unique(unlist(strsplit(bad, "")))
  code <- paste("%", charToRaw(paste(bad, collapse="")), sep="")
  for (i in seq_along(bad))
    str <- gsub(bad[i], code[i], str, fixed = TRUE)
  str
}

# differs from URLdecode (vectorized)
urlDecode <- function(str)
{
  unistr <- paste(str, collapse="")
  percents <- gregexpr("%", unistr, fixed=TRUE)
  pos <- percents[[1]]
  if (pos[1] != -1) {
    code <- unique(substring(unistr, pos, pos+2))
    raws <- as.raw(gsub("%", "0x", code, fixed=TRUE))
    real <- strsplit(rawToChar(raws), "")[[1]]
    for (i in seq_along(code))
      str <- gsub(code[i], real[i], str, fixed=TRUE)
  }
  str
}
