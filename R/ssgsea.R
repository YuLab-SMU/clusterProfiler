# This function performs ssGSEA Analysis



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

# ssgsea(gene_expression, genes, num_gene_sets, genes_per_set)
