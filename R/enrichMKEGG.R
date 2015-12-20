##' KEGG Module Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment KEGG Module
##' categories with FDR control.
##'
##'
##' @inheritParams enrichKEGG
##' @return A \code{enrichResult} instance.
##' @export
enrichMKEGG <- function(gene,
                        species = 'hsa',
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        universe,
                        minGSSize = 5,
                        qvalueCutoff = 0.2) {
    
    KEGG_DATA <- download.KEGG(species, "MKEGG")
    res <- enricher_internal(gene,
                             pvalueCutoff  =pvalueCutoff,
                             pAdjustMethod =pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             USER_DATA = KEGG_DATA)
    
    res@ontology <- "MKEGG"
    res@organism <- species
    return(res)
}
