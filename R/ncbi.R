NCBI_TAX_URL <- "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"

isNCBISpeciesURL <- function(url) {
    substring(url, 1L, nchar(NCBI_TAX_URL)) == NCBI_TAX_URL
}

metadataFromNCBI <- function(url) {
    html <- httpGet(url)
    list("Taxonomy ID"=parseTaxonomyIdFromNCBI(html),
         "Organism"=parseOrganismFromNCBI(html))
}

parseTaxonomyIdFromNCBI <- function(html) {
    unlist(getNodeSet(html, "//input[@name='old_id']/@value"), use.names=FALSE)
}

taxonomyIdFromNCBI <- function(species) {
    url <- paste0(NCBI_TAX_URL, "?name=", URLencode(species))
    parseTaxonomyIdFromNCBI(httpGet(url))
}

parseOrganismFromNCBI <- function(html) {
    title <- xmlValue(getNodeSet(html, "//title/text()")[[1L]])
    sub(".*\\((.*?)\\).*", "\\1", title)
}

speciesFromNCBI <- function(id) {
    url <- paste0(NCBI_TAX_URL, "?id=", URLencode(id))
    parseOrganismFromNCBI(httpGet(url))
}
