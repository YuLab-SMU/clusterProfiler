#' A universal enrichment analyzer
#'
#'
#' @title enricher
#' @param gene a vector of gene id
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
#' @param minGSSize minimal size of genes annotated for testing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant.  Tests must pass i) \code{pvalueCutoff} on unadjusted pvalues, ii) \code{pvalueCutoff} on adjusted pvalues and iii) \code{qvalueCutoff} on qvalues to be reported.
#' @param gson a GSON object, if not NULL, use it as annotation data. 
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene.
#' Only used when gson is NULL.
#' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name.
#' Only used when gson is NULL.
#' @return A \code{enrichResult} instance
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @export
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
                enrichit::ora_gson(gene = gene,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      universe = universe,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      qvalueCutoff = qvalueCutoff,
                      gson = USER_DATA)
        })
        class(res) <- "enrichResultList"
        return(res)
    }

    if (is.null(gson)) {
        # USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
        gsid2gene <- TERM2GENE
        colnames(gsid2gene) <- c("gsid", "gene")
        
        if (missing(TERM2NAME) || is.null(TERM2NAME) || all(is.na(TERM2NAME))) {
            gsid2name <- data.frame(gsid=unique(gsid2gene$gsid), name=unique(gsid2gene$gsid))
        } else {
            gsid2name <- TERM2NAME
            colnames(gsid2name) <- c("gsid", "name")
        }
        
        USER_DATA <- gson::gson(gsid2gene = gsid2gene,
                                gsid2name = gsid2name,
                                species = "unknown",
                                gsname = "unknown",
                                version = "unknown",
                                accessed_date = as.character(Sys.Date()),
                                keytype = "unknown")
    } else {
        if (!inherits(gson,  "GSON")) {
            stop("gson shoud be a GSON or GSONList object")
        }
        USER_DATA <- gson
    }
    
    enrichit::ora_gson(gene = gene,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      universe = universe,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      qvalueCutoff = qvalueCutoff,
                      gson = USER_DATA)
}


#' a universal gene set enrichment analysis tools
#'
#'
#' @title GSEA
#' @param geneList order ranked geneList
#' @param exponent weight of each step
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param nPerm The number of permutations.
#' @param method method of calculating the pvalue, one of "multilevel", "monte carlo" and "fgsea"
#' @param adaptive logical, whether to use adaptive method for calculating pvalue
#' @param minPerm minimal number of permutations for adaptive method
#' @param maxPerm maximal number of permutations for adaptive method
#' @param pvalThreshold pvalue threshold for adaptive method
#' @param eps boundary for calculating the p value in multilevel mode
#' @param pvalueCutoff p-value cutoff applied to both the raw p-value and the adjusted p-value (p.adjust), consistent with the historical clusterProfiler/DOSE behavior
#' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param gson a GSON object, if not NULL, use it as annotation data. 
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene.
#' Only used when gson is NULL.
#' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name.
#' Only used when gson is NULL.
#' @param verbose logical
#' @param seed random seed for reproducibility, set to a number (or TRUE to use a
#'   fixed default seed) to make the result reproducible, or FALSE (default) to draw
#'   a random seed on each run, so results may vary between runs. The underlying
#'   permutation engine uses its own RNG seeded with this value; see
#'   \code{enrichit::gsea()} for details.
#' @param ... other parameters passed to \code{enrichit::gsea_gson()}
#' @return gseaResult object
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @export
GSEA <- function(geneList,
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 1e-10,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 gson  = NULL,
                 TERM2GENE,
                 TERM2NAME = NA,
                 verbose = TRUE,
                 nPerm = 1000,
                 method = "multilevel",
                 adaptive = FALSE,
                 minPerm = 101,
                 maxPerm = 1e5,
                 pvalThreshold = 0.1,
                 seed = FALSE,
                 ...) {

    if (inherits(gson, 'GSONList')) {
        res <- lapply(gson, function(USER_DATA) {
            enrichit::gsea_gson(geneList      = geneList,
                          exponent      = exponent,
                          minGSSize     = minGSSize,
                          maxGSSize     = maxGSSize,
                          eps           = eps,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose       = verbose,
                          gson          = USER_DATA,
                          nPerm         = nPerm,
                          method        = method,
                          adaptive      = adaptive,
                          minPerm       = minPerm,
                          maxPerm       = maxPerm,
                          pvalThreshold = pvalThreshold,
                          seed          = seed,
                          ...)
        })
        
        class(res) <- "gseaResultList"
        return(res)
    }
    if (is.null(gson)) {
        # USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
        gsid2gene <- TERM2GENE
        colnames(gsid2gene) <- c("gsid", "gene")
        
        if (missing(TERM2NAME) || is.null(TERM2NAME) || all(is.na(TERM2NAME))) {
            gsid2name <- data.frame(gsid=unique(gsid2gene$gsid), name=unique(gsid2gene$gsid))
        } else {
            gsid2name <- TERM2NAME
            colnames(gsid2name) <- c("gsid", "name")
        }
        
        USER_DATA <- gson::gson(gsid2gene = gsid2gene,
                                gsid2name = gsid2name,
                                species = "unknown",
                                gsname = "unknown",
                                version = "unknown",
                                accessed_date = as.character(Sys.Date()),
                                keytype = "unknown")
    } else {
        if (!inherits(gson,  "GSON")) {
            stop("gson shoud be a GSON or GSONList object")
        }
        USER_DATA <- gson
    }
    
    enrichit::gsea_gson(geneList      = geneList,
                  exponent      = exponent,
                  minGSSize     = minGSSize,
                  maxGSSize     = maxGSSize,
                  eps           = eps,
                  pvalueCutoff  = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  verbose       = verbose,
                  gson          = USER_DATA,
                  nPerm         = nPerm,
                  method        = method,
                  adaptive      = adaptive,
                  minPerm       = minPerm,
                  maxPerm       = maxPerm,
                  pvalThreshold = pvalThreshold,
                  seed          = seed,
                  ...)
    
}

