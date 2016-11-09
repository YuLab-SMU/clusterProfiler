########### zzz.R
.onLoad <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields="Version")
    msg <- paste0(pkgname, " v", pkgVersion, "\n",
                  "For help: https://guangchuangyu.github.io/clusterProfiler")
    packageStartupMessage(msg)
    .initial()

}


.initial <- function() {
    ## assign(".clusterProfilesEnv", new.env(),.GlobalEnv)
}
