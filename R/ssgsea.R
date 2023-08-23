##' Single Sample Gene Set Enrichment Analysis (ssGSEA)
##'
##' This function performs single-sample gene set enrichment analysis (ssGSEA) 
##'
##' @title ssgsea
##' @param gene_expression a matrix of gene expression data with genes in rows and samples in columns
##' @param gene_sets a list of gene sets used for enrichment analysis.
##' @param num_gene_sets the number of gene sets to generate and use for ssGSEA (default: 100).
##' @return A matrix of ssGSEA enrichment scores with gene sets as rows and samples as columns.
##' @export

ssgsea<- function(gene_expression, gene_sets, num_gene_sets=100) {
  gene_expression <- convert_id(gene_expression)
  gene_sets <- gene_sets[1:num_gene_sets]
  ssGSEA_scores <- gsva(gene_expression,gene_sets, method = "ssgsea")
  return(ssGSEA_scores)
}

convert_id <- function(gene_expr) {
  genes <- rownames(gene_expr)
  entrez_ids <- mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
  uniprot_ids <- mapIds(org.Hs.eg.db, keys=entrez_ids, column="UNIPROT", keytype="ENTREZID", multiVals="first")
  rownames(gene_expr) <- uniprot_ids
  converted_genes <- rownames(gene_expr)
  return(gene_expr)
}


