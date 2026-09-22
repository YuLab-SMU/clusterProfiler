#' KEGG Module Enrichment Analysis of a gene set.
#' Given a vector of genes, this function will return the enrichment KEGG Module
#' categories with FDR control.
#'
#' As with \code{enrichKEGG()}, \code{organism} may be a species code or a
#' \code{GSON} object, so module enrichment can be run against a locally built
#' annotation (see \code{gson_KEGG()}) when KEGG is unreachable.
#'
#' @inheritParams enrichKEGG
#' @return A \code{enrichResult} instance.
#' @export
enrichMKEGG <- function(gene,
                        organism = 'hsa',
                        keyType = 'kegg',
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        universe,
                        minGSSize = 10,
                        maxGSSize = 500,
                        qvalueCutoff = 0.2) {

    if (inherits(organism, "GSON")) {
        KEGG_DATA <- organism
        species <- KEGG_DATA@species
        keyType <- KEGG_DATA@keytype
    } else if (inherits(organism, "character")) {
        species <- organismMapper(organism)
        KEGG_DATA <- prepare_KEGG(species, "MKEGG", keyType)
    } else {
        stop("organism should be a species name or a GSON object")
    }

    res <- enrichit::ora_gson(gene,
                             pvalueCutoff  = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             maxGSSize     = maxGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             gson = KEGG_DATA)

    if (is.null(res))
        return(res)
    
    
    res@ontology <- "MKEGG"
    res@organism <- species
    ## a GSON knows its own keytype; a species name does not, and reporting
    ## "UNKNOWN" there is what this function has always done
    res@keytype <- if (inherits(organism, "GSON")) keyType else "UNKNOWN"
    
    return(res)
}
