
.ConnectionManager <- setRefClass("ConnectionManager",
                                  fields = c(connections = "list"))

manager <- function() .ConnectionManager()

connection <- function(manager, x, open = "") {
  connectionForResource(manager, resource(x), open = open)
}

resourceDescription <- function(x) {
  r <- resource(x)
  if (is(r, "connection"))
    r <- summary(r)$description
  r
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

## First checks for Windows drive letter.
## There are no known URI schemes that are only a single character.
isURL <- function(uri) {
    if (!isSingleString(uri))
        return(FALSE)
    windowsDriveLetter <- .Platform$OS.type == "windows" &&
        grepl("^[A-Za-z]:[/\\]", uri)
    grepl("^[A-Za-z]+:", uri) && !windowsDriveLetter
}

## Uses XML::parseURI, except custom check for whether it is a URL
.parseURI <- function(uri) {
  if (!isURL(uri)) {
    parsed <- parseURI("")
    parsed$path <- uri
  } else {
    parsed <- parseURI(uri)
    if (parsed$scheme == "file" && .Platform$OS.type == "windows") 
      parsed$path <- substring(parsed$path, 2) # trim '/' from '/C:/foo/bar.txt'
  }
  parsed
}

normURI <- function(x) {
  if (!isSingleString(x))
    stop("URI must be a single, non-NA string")
  uri <- .parseURI(x)
  if (uri$scheme == "") # /// (vs. //) needed for Windows
    x <- paste("file:///", file_path_as_absolute(x), sep = "")
  x
}

createResource <- function(x, dir = FALSE, content = "") {
  uri <- .parseURI(x)
  if (uri$scheme == "file" || uri$scheme == "") {
    if (!file.exists(uri$path)) {
      if (dir)
        dir.create(uri$path, recursive = TRUE)
      else writeLines(content, uri$path)
    } else warning("Path '", uri$path, "' already exists")
  } else stop("Cannot create a resource that is not a local file")
}

uriExists <- function(x) {
  uri <- .parseURI(x)
  if (uriIsLocal(x)) {
    exists <- file.exists(uri$path)
  } else {
    txt <- getURL(x, header = TRUE)
    exists <- grepl("^HTTP/\\d+\\.\\d+ 200 OK", txt)
  }
  exists
}

uriIsLocal <- function(x) {
  x$scheme == "file" || x$scheme == ""
}

uriIsWritable <- function(x) {
  uri <- .parseURI(x)
  if (uriIsLocal(uri)) {
    !file.access(uri$path, 2) ||
    (!file.exists(uri$path) && uriIsWritable(dirname(uri$path)))
  } else FALSE
}

checkArgFormat <- function(con, format) {
  if (toupper(format) !=
      substring(toupper(sub("File$", "", class(con))), 1, nchar(format)))
    stop("Cannot treat a '", class(con), "' as format '", format, "'")
}

connectionForResource <- function(manager, x, open = "") {
  resource <- decompress(manager, x)
  if (is.character(resource)) {
    if (!nzchar(resource))
      stop("path cannot be an empty string")
    uri <- .parseURI(resource)
    if (uri$scheme != "")
      con <- url(resource)
    else con <- file(resource)
  } else con <- resource
  if (!isOpen(con) && nzchar(open)) {
      open(con, open)
      con <- manage(manager, con)
  }
  con
}

## Connection management (similar to memory management)

manage <- function(manager, con) {
    manager$connections <- unique(c(manager$connections, list(con)))
    attr(con, "manager") <- manager
    con
}

managed <- function(manager, con) {
    con %in% manager$connections
}

unmanage <- function(manager, con) {
    manager$connections <- setdiff(manager$connections, con)
    attr(con, "manager") <- NULL
    con
}

release <- function(manager, con) {
    if (managed(manager, con)) {
        unmanage(manager, con)
        close(con)
    }
    con
}
