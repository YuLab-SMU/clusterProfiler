##' Single Sample Gene Set Enrichment Analysis (ssGSEA)
##'
##' This function performs single-sample gene set enrichment analysis (ssGSEA) 
##'
##' @title ssgsea
##' @param gene_expression a matrix of gene expression data with genes in rows and samples in columns
##' @param genes a vector of gene identifiers
##' @param num_gene_sets the number of gene sets to generate and use for ssGSEA (default: 100)
##' @param genes_per_set the number of genes to include in each generated gene set (default: 10)
##' @return A matrix of ssGSEA enrichment scores with gene sets as rows and samples as columns
##' @export

ssgsea <- function(gene_expression, genes, num_gene_sets=100, genes_per_set=10) {
  gs <- list()
  for (i in 1:num_gene_sets) {
    sampled_genes <- sample(genes, genes_per_set, replace = FALSE)
    gs[[i]] <- sampled_genes
  }
  names(gs) <- paste0("gs", 1:length(gs))

  ssGSEA_scores <- gsva(gene_expression, gs, method = "ssgsea")
  return(ssGSEA_scores)
}
