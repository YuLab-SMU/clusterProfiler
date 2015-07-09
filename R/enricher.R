##' A universal enrichment analyzer 
##'
##' 
##' @title enricher
##' @param gene a vector of gene id
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff
##' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
##' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name
##' @return A \code{enrichResult} instance
##' @author Guangchuang Yu
##' @export
enricher <- function(gene,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe,
                     minGSSize=5,
                     qvalueCutoff = 0.2,
                     TERM2GENE,
                     TERM2NAME = NA) {
    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
    enrich.internal(gene = gene,
                    organism = "UNKNOWN",
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod = pAdjustMethod,
                    ont = "USER_DEFINED",
                    universe = universe,
                    minGSSize = minGSSize,
                    qvalueCutoff = qvalueCutoff,
                    readable = FALSE,
                    USER_DATA = USER_DATA)
}
                     
                     
##' a universal gene set enrichment analysis tools
##'
##' 
##' @title GSEA
##' @param geneList order ranked geneList 
##' @param exponent weight of each step
##' @param nPerm number of permutations
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod p value adjustment method
##' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
##' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name
##' @param verbose logical
##' @return gseaResult object
##' @author Guangchuang Yu
##' @export
GSEA <- function(geneList,
                 exponent = 1,
                 nPerm = 1000,
                 minGSSize = 10,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 TERM2GENE,
                 TERM2NAME = NA,
                 verbose = TRUE) {
    if(verbose)
        cat("preparing geneSet collections...\n")

    organism = "UNKNOWN"
    setType <- "USER_DEFINED"
    class(setType) <- setType
    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
    geneSets <- getGeneSet(setType, organism, USER_DATA = USER_DATA)

    gsea(geneList = geneList,
         geneSets = geneSets,
         setType = setType,
         organism = organism,
         exponent = exponent,
         nPerm = nPerm,
         minGSSize = minGSSize,
         pvalueCutoff = pvalueCutoff,
         pAdjustMethod = pAdjustMethod,
         verbose = verbose,
         USER_DATA = USER_DATA)
}

##' @method getGeneSet USER_DEFINED
##' @export
getGeneSet.USER_DEFINED <- function(setType = "USER_DEFINED", organism, ...) {
    getGeneSet.USER_DEFINED.internal(setType, organism, ...)
}

getGeneSet.USER_DEFINED.internal <- function(setType, organism, USER_DATA, ...) {
    gs <- get("PATHID2EXTID", envir = USER_DATA)
    return(gs)
}

