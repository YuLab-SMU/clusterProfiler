##' parse gmt file to a data.frame
##'
##'
##' @title read.gmt
##' @param gmtfile gmt file
## @importFrom GSEABase getGmt
## @importFrom GSEABase geneIds
##' @importFrom utils stack
##' @importFrom rvcheck get_fun_from_pkg
##' @export
##' @return data.frame
##' @author Guangchuang Yu
read.gmt <- function(gmtfile) {
    getGmt <- get_fun_from_pkg("GSEABase", "getGmt")
    geneIds <- get_fun_from_pkg("GSEABase", "geneIds")

    gmt <- getGmt(con=gmtfile)
    ont2gene <- geneIds(gmt) %>% stack
    ont2gene <- ont2gene[, c("ind", "values")]
    colnames(ont2gene) <- c("ont", "gene")
    return(ont2gene)
}

