##' A universal enrichment analyzer
##'
##'
##' @title enricher
##' @param gene a vector of gene id
##' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
##' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
##' @param minGSSize minimal size of genes annotated for testing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant.  Tests must pass i) \code{pvalueCutoff} on unadjusted pvalues, ii) \code{pvalueCutoff} on adjusted pvalues and iii) \code{qvalueCutoff} on qvalues to be reported.
##' @param gson a GSON object, if not NULL, use it as annotation data. 
##' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene.
##' Only used when gson is NULL.
##' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name.
##' Only used when gson is NULL.
##' @return A \code{enrichResult} instance
##' @author Guangchuang Yu \url{https://yulab-smu.top}
##' @export
enricher <- function(gene,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = NULL,
                     minGSSize=10,
                     maxGSSize=500,
                     qvalueCutoff = 0.2,
                     gson  = NULL,
                     TERM2GENE,
                     TERM2NAME = NA
                     ) {
    if (inherits(gson, 'GSONList')) {
        res <- lapply(gson, function(USER_DATA) {
                enricher_internal(gene = gene,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      universe = universe,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      qvalueCutoff = qvalueCutoff,
                      USER_DATA = USER_DATA)
        })
        class(res) <- "enrichResultList"
        return(res)
    }

    if (is.null(gson)) {
        USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
    } else {
        if (!inherits(gson,  "GSON")) {
            stop("gson shoud be a GSON or GSONList object")
        }
        USER_DATA <- gson
    }
    
    enricher_internal(gene = gene,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      universe = universe,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      qvalueCutoff = qvalueCutoff,
                      USER_DATA = USER_DATA)
}


##' a universal gene set enrichment analysis tools
##'
##'
##' @title GSEA
##' @param geneList order ranked geneList
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param eps This parameter sets the boundary for calculating the p value.
##' @param pvalueCutoff adjusted pvalue cutoff
##' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param gson a GSON object, if not NULL, use it as annotation data. 
##' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene.
##' Only used when gson is NULL.
##' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name.
##' Only used when gson is NULL.
##' @param verbose logical
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @param ... other parameter
##' @return gseaResult object
##' @author Guangchuang Yu \url{https://yulab-smu.top}
##' @export
GSEA <- function(geneList,
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps  = 1e-10,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 gson  = NULL,
                 TERM2GENE,
                 TERM2NAME = NA,
                 verbose = TRUE,
                 seed = FALSE,
                 by = 'fgsea',
                 ...) {

    if (inherits(gson, 'GSONList')) {
        res <- lapply(gson, function(USER_DATA) {
            GSEA_internal(geneList      = geneList,
                          exponent      = exponent,
                          minGSSize     = minGSSize,
                          maxGSSize     = maxGSSize,
                          eps           = eps,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose       = verbose,
                          USER_DATA     = USER_DATA,
                          seed          = seed,
                          by            = by,
                          ...)
        })
        
        class(res) <- "gseaResultList"
        return(res)
    }
    if (is.null(gson)) {
        USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
    } else {
        if (!inherits(gson,  "GSON")) {
            stop("gson shoud be a GSON or GSONList object")
        }
        USER_DATA <- gson
    }
    
    GSEA_internal(geneList      = geneList,
                  exponent      = exponent,
                  minGSSize     = minGSSize,
                  maxGSSize     = maxGSSize,
                  eps           = eps,
                  pvalueCutoff  = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  verbose       = verbose,
                  USER_DATA     = USER_DATA,
                  seed          = seed,
                  by            = by,
                  ...)
    
}

