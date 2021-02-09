#' get universe microbiome data
#'
#' @param Kgenelist a list of K gene numbers of your microbiome data as universe
#'
#' @return A gene annotation table of  human gut microbiome
#' @export
#'
#' @examples
#' \dontrun{
#' getM_DATA(hgmlist)
#' }
#'
#'
getM_DATA <- function(Kgenelist){
    p <- kegg_list('pathway')
    p2 <- kegg_link('ko',"pathway")
    p2 <- p2[grep(pattern="path:map",p2[,1]),]
    res <- merge(p2, p, by = 'from', all.x=TRUE)
    colnames(res) <- c("pathway", "Knum", "name")
    res$Knum <- gsub("ko:","",res$Knum)
    M_DATA <- res[res[,2] %in% Kgenelist,]
    return(M_DATA)
}



#' KEGG enrichment analysis for microbiome.
#'
#' @param gene a vector of K gene id (e.g. K00001).
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param universe universe background genes. If missing, the all genes listed in the microbe_data will be used as background.
#' @param minGSSize minGSSize minimal size of genes annotated by KEGG term for testing.
#' @param maxGSSize maxGSSize maximal size of genes annotated for testing.
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant.
#' @param use_internal_data
#' @param microbe_data a list of all K gene numbers of your microbiome data as universe.
#'
#' @return A \code{enrichResult} instance.
#' @export
#'
#' @examples
#' \dontrun{
#' data(hgmlist)
#' yy <- enrichmbKEGG(hgmlist[1:20],microbe_data=hgmlist)
#' }
#'
#'
enrichmbKEGG <- function(gene,
                         pvalueCutoff      = 0.05,
                         pAdjustMethod     = "BH",
                         universe,
                         minGSSize         = 10,
                         maxGSSize         = 500,
                         qvalueCutoff      = 0.2,
                         use_internal_data = FALSE,
                         microbe_data) {
    M_DATA <- getM_DATA(microbe_data)
    KEGG_DATA <- build_Anno(M_DATA[c(c("pathway","Knum"))], M_DATA[c("pathway","name")])
    res <- enricher_internal(gene,
                             pvalueCutoff  = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             maxGSSize     = maxGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             USER_DATA = KEGG_DATA)
    if (is.null(res))
        return(res)

    res@ontology <- "KEGG"
    res@organism <- "microbiome"

    return(res)
}



