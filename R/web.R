htmlErrorHandler <- function(msg, code, domain, line, col, level, filename) {
  if (!length(level))
    stop("Unknown HTML parse error")
  if (level > 2)
    stop("Failed to Parse HTML [", line, ":", col, "]: ", msg)
}

htmlParse <- function(str)
  suppressWarnings(htmlTreeParse(str, asText = TRUE, useInternalNodes = TRUE,
                                 error = htmlErrorHandler))

httpGet <- function(url, .form = list(), .parse = TRUE, ...) {
  if (length(.form) == 0)
      out <- getURL(url, useragent = "rtracklayer",
                    verbose=getOption("rtracklayer.http.verbose", FALSE), ...)
  else out <- getForm(url, .params = .form,
                      .opts=list(..., useragent = "rtracklayer",
                                 verbose=getOption("rtracklayer.http.verbose",
                                                   FALSE)))
  if (.parse)
    htmlParse(out)
  else out
}

httpPost <- function(url, .form = list(), .parse = TRUE, ...) {
  form <- postForm(url, .params = .form,
                   .opts=list(..., useragent = "rtracklayer",
                              verbose=getOption("rtracklayer.http.verbose",
                                                FALSE)))
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
  str <- as.character(str)
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
urlDecode <- function(str, na.strings="NA")
{
  ans <- curlUnescape(str)
  if (!identical(na.strings, "NA"))
      ans[is.na(str)] <- na.strings
  ans
}

expandPath <- function(x) {
    if (startsWith(x, "http") || startsWith(x, "ftp"))
        expandURL(x)
    else path.expand(x)
}

expandURL <- function(uri) {
    if(!url.exists(uri))
        return(uri)
    else {
        opts <- list(
            followlocation = TRUE,  # resolve redirects
            ssl.verifyhost = FALSE, # suppress certain SSL errors
            ssl.verifypeer = FALSE, 
            nobody = TRUE, # perform HEAD request
            verbose = FALSE
            )
        curlhandle <- getCurlHandle(.opts = opts)
        getURL(uri, curl = curlhandle)
        info <- getCurlInfo(curlhandle)
        info$effective.url
    }
}
