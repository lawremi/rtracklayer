.onUnload <- function(libpath)
{
    library.dynam.unload("rtracklayer", libpath)
}

setUserUdcDir <- function() {
    dir <- paste0("/tmp/udcCache_", Sys.info()[["user"]])
    .Call(R_setUserUdcDir, dir)
}

.onLoad <- function(libname, pkgname)
{
    setUserUdcDir()
    S7::methods_register()
    ## S7/issues/540
    rm(list = c("[[", "[[<-", "names"), envir = getNamespace("rtracklayer"))
}
