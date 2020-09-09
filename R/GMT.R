##' parse gmt file to a data.frame
##'
##'
##' @title read.gmt
##' @rdname read-gmt
##' @param gmtfile gmt file
##' @importFrom utils stack
##' @importFrom rvcheck get_fun_from_pkg
##' @export
##' @return data.frame
##' @author Guangchuang Yu
read.gmt <- function(gmtfile) {
    ## getGmt <- get_fun_from_pkg("GSEABase", "getGmt")
    ## geneIds <- get_fun_from_pkg("GSEABase", "geneIds")

    ## gmt <- getGmt(con=gmtfile)
    ## ont2gene <- geneIds(gmt) %>% stack
    ## ont2gene <- ont2gene[, c("ind", "values")]
    ## colnames(ont2gene) <- c("ont", "gene")

    x <- readLines(gmtfile)
    ## ont2gene <- lapply(x, function(record) {
    ##     y = strsplit(record, "\t")[[1]]
    ##     data.frame(ont=y[1], gene=y[-c(1:2)])
    ## }) %>% do.call('rbind', .)
    

    ## first column: gene set name
    ## second column: description
    ## all the others, unequal length for genes
    res <- strsplit(x, "\t")
    names(res) <- vapply(res, function(y) y[1], character(1))
    res <- lapply(res, "[", -c(1:2))

    ont2gene <- stack(res)
    ont2gene <- ont2gene[, c("ind", "values")]
    colnames(ont2gene) <- c("term", "gene")
    return(ont2gene)
}

##' @rdname read-gmt
##' @importFrom rlang .data
read.gmt.wp <- function(gmtfile) {
    read.gmt(gmtfile) %>%
        tidyr::separate(.data$term, c("name","version","wpid","org"), "%")
}

