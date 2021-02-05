##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields="Version")
    msg <- paste0(pkgname, " v", pkgVersion, "  ",
                  "For help: https://guangchuangyu.github.io/software/", pkgname, "\n\n")

    if (capabilities("libcurl")) {
        dl.method <- "libcurl"
    } else {
        dl.method <- getOption("download.file.method", default = "auto")        
    }

    options(clusterProfiler.download.method = dl.method)
    options(timeout = max(300, getOption("timeout"))) # see ?download.file
    
    citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                       "Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. ",
                       "clusterProfiler: an R package for comparing biological themes among gene clusters. ",
                       "OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.")

    packageStartupMessage(paste0(msg, citation))
}


