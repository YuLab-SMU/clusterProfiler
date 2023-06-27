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
    res@keytype <-  "HGNCID"

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
gsePC <- function(geneList, organism, ...) {
    pcdata <- prepare_PC_data(organism)
    res <- GSEA(geneList,
                TERM2GENE = pcdata$PCID2GENE,
                TERM2NAME = pcdata$PCID2NAME,
                ...)

    if (is.null(res)) return(res)

    res@setType <- "Pathway Commons"
    res@organism <- organism
    res@keytype <-  "HGNCID"

    return(res)
}

##' @importFrom rlang .data
prepare_PC_data <- function(organism) {
  pc2gene <- get_pc_data(organism)
  ##TERM2GENE
  pcid2gene <- pc2gene %>% dplyr::select(.data$pcid, .data$gene)
  ##TERM2NAME
  pcid2name <- pc2gene %>% dplyr::select(.data$pcid, .data$name)
  list(PCID2GENE = pcid2gene,
       PCID2NAME = pcid2name)
}

get_pc_gmtfile <- function() {
  pcurl <- 'https://www.pathwaycommons.org/archives/PC2/v12/'
  x <- readLines(pcurl)
  y <- x[grep('\\.gmt',x)]
  sub(".*(PathwayCommons.*\\.gmt.gz).*", "\\1",  y[grep('', y)])
}

get_pc_source <- function() {
    gmtfile <- get_wp_gmtfile()
    orgs <- sub("PathwayCommons\\d+\\.([_A-Za-z]+)\\.([_A-Za-z]+)\\.gmt.gz", "\\1", gmtfile)
}

