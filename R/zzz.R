##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields="Version")
    msg <- paste0(pkgname, " v", pkgVersion, "  ",
                  "For help: https://yulab-smu.top/biomedical-knowledge-mining-book/", "\n\n")
    if (capabilities("libcurl")) {
        dl.method <- "libcurl"
    } else {
        dl.method <- getOption("download.file.method", default = "auto")        
    }

    options(clusterProfiler.download.method = dl.method)
    options(timeout = max(300, getOption("timeout"))) # see ?download.file
    
    citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                       "Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L. ", 
                       "Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., Yu, G. ", 
                       "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ", 
                       "The Innovation. 2021, doi: 10.1016/j.xinn.2021.100141.")
    packageStartupMessage(paste0(msg, citation))
}


