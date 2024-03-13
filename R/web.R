htmlErrorHandler <- function(msg, code, domain, line, col, level, filename) {
  if (!length(level))
    stop("Unknown HTML parse error")
  if (level > 2)
    stop("Failed to Parse HTML [", line, ":", col, "]: ", msg)
}

htmlParse <- function(str)
  suppressWarnings(htmlTreeParse(str, asText = TRUE, useInternalNodes = TRUE,
                                 error = htmlErrorHandler))

rtracklayerGET <- function(url, ..., query = list()) {
  verbose <- getOption("rtracklayer.http.verbose", FALSE)
  verbose <- as.integer(isTRUE(as.logical(verbose)))
  response <- GET(url, user_agent("rtracklayer"), config(verbose = verbose),
                       ...,
                       query = as.list(query))
  htmlParse(content(response, as="text"))
}

rtracklayerPOST <- function(url, ..., body = list()) {
  verbose <- getOption("rtracklayer.http.verbose", FALSE)
  verbose <- as.integer(isTRUE(as.logical(verbose)))
  response <- POST(url, user_agent("rtracklayer"), config(verbose = verbose),
                        ...,
                        body = as.list(body))
  htmlParse(content(response, as = "text"))
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

urlDecode <- function(str, na.strings="NA")
{
  ans <- URLdecode(str)
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
    if(HEAD(uri)$status_code != 200L)
        return(uri)
    else {
        #opts <- list(
        #    followlocation = TRUE,  # resolve redirects
        #    ssl.verifyhost = FALSE, # suppress certain SSL errors
        #    ssl.verifypeer = FALSE,
        #    nobody = TRUE, # perform HEAD request
        #    verbose = FALSE
        #    )
        #curlhandle <- getCurlHandle(.opts = opts)
        #getURL(uri, curl = curlhandle)
        #info <- getCurlInfo(curlhandle)
        #info$effective.url

        config <- list(followlocation = 1,
                       ssl_verifyhost = 0,
                       ssl_verifypeer = 0,
                       nobody = 1,
                       verbose = 0)
        ## nobody = 1 seems to be ignored and HEAD() doesn't work here
        GET(uri, config)$url
    }
}
