##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary

##' @importFrom DOSE geneID
##' @export
DOSE::geneID

##' @importFrom DOSE geneInCategory
##' @export
DOSE::geneInCategory

##' @importFrom DOSE gsfilter
##' @export
DOSE::gsfilter

##' @importFrom DOSE setReadable
##' @export
DOSE::setReadable

build_Anno <- getFromNamespace("build_Anno", "DOSE")
TERMID2EXTID <- getFromNamespace("TERMID2EXTID", "DOSE")
TERM2NAME <- getFromNamespace("TERM2NAME", "DOSE")
get_organism <- getFromNamespace("get_organism", "DOSE")
enricher_internal <- getFromNamespace("enricher_internal", "DOSE")

