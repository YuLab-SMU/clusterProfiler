#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    # pkgVersion <- packageDescription(pkgname, fields="Version")
    # msg <- paste0(pkgname, " v", pkgVersion, "  ",
    #               "For help: https://yulab-smu.top/biomedical-knowledge-mining-book/", "\n\n")
    
    # A per-session cache used to default to tempdir(), unconditionally. That
    # overwrote any directory the user set in .Rprofile, and because tempdir()
    # differs in every session it also made downloaded resources -- the
    # mechanism_sets() GMT files -- miss their cache and download again each
    # time. Only fill the option when nobody has set it.
    if (is.null(getOption("clusterProfiler_cache_dir"))) {
        options(clusterProfiler_cache_dir = yulab.utils::user_dir("clusterProfiler"))
    }

    packageStartupMessage(yulab.utils::yulab_msg(pkgname))
}


