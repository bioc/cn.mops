# TODO: Add comment
# 
# Author: klambaue
###############################################################################


.onLoad <- function(libname, pkgname) {
	library.dynam("cn.mops", pkgname, libname)
}

.onUnload <- function(libpath)
{
	library.dynam.unload("cn.mops", libpath)
}

.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.12")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
