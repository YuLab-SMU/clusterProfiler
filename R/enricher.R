##' A universal enrichment analyzer
##'
##'
##' @title enricher
##' @param gene a vector of gene id
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated for testing
##' @param maxGSSize maximal size of genes annotated for testing
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
                     minGSSize=10,
                     maxGSSize=500,
                     qvalueCutoff = 0.2,
                     TERM2GENE,
                     TERM2NAME = NA) {
    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
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
##' @param nPerm number of permutations
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod p value adjustment method
##' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
##' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name
##' @param verbose logical
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @return gseaResult object
##' @author Guangchuang Yu
##' @export
GSEA <- function(geneList,
                 exponent = 1,
                 nPerm = 1000,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 TERM2GENE,
                 TERM2NAME = NA,
                 verbose = TRUE,
                 seed = FALSE,
                 by = 'fgsea') {

    USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)

    GSEA_internal(geneList = geneList,
                  exponent = exponent,
                  nPerm = nPerm,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  verbose = verbose,
                  USER_DATA = USER_DATA,
                  seed = seed,
                  by = by)
}

