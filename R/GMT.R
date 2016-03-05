##' parse gmt file to a data.frame
##'
##' 
##' @title read.gmt
##' @param gmtfile gmt file
##' @importFrom GSEABase getGmt
##' @importFrom GSEABase geneIds
##' @importFrom utils stack
##' @export
##' @return data.frame 
##' @author Guangchuang Yu
read.gmt <- function(gmtfile) {
    gmt <- getGmt(con=gmtfile)
    ont2gene <- geneIds(gmt) %>% stack
    ont2gene <- ont2gene[, c("ind", "values")]
    colnames(ont2gene) <- c("ont", "gene")
    return(ont2gene)
}

