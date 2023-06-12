##' ORA analysis for Pathway Commons
##'
##' This function performs over-representation analysis using  Pathway Commons
##' @title enrichPC
##' @param gene a vector of entrez gene id
##' @param organism supported organisms, which can be accessed via the get_pc_organisms() function
##' @param ... additional parameters, see also the parameters supported by the enricher() function
##' @return A \code{enrichResult} instance
##' @export

enrichPC <- function(gene, organism, ...) {
    pcdata <- prepare_PC_data(organism)
    res <- enricher(gene,
                    TERM2GENE = pcdata$PCID2GENE,
                    TERM2NAME = pcdata$PCID2NAME,
                    ...)
    if (is.null(res)) return(res)

    res@ontology <- "Pathway Commons"
    res@organism <- organism
    res@keytype <-  "ENTREZID"

    return(res)
}

##' GSEA analysis for  Pathway Commons
##'
##' This function performs GSEA using  Pathway Commons
##' @title gsePC
##' @param geneList ranked gene list
##' @param organism supported organisms, which can be accessed via the get_pc_organisms() function
##' @param ... additional parameters, see also the parameters supported by the GSEA() function
##' @return A \code{gseaResult} instance
##' @export
##' @author Guangchuang Yu 
gsePC <- function(geneList, organism, ...) {
    pcdata <- prepare_PC_data(organism)
    res <- GSEA(geneList,
                TERM2GENE = pcdata$PCID2GENE,
                TERM2NAME = pcdata$PCID2NAME,
                ...)

    if (is.null(res)) return(res)

    res@setType <- "Pathway Commons"
    res@organism <- organism
    res@keytype <-  "ENTREZID"

    return(res)
}

