library(testthat)

make_kegg_enrich <- function() {
  df <- data.frame(
    ID = c("hsa04510", "hsa04512", "hsa04110", "hsa03030"),
    Description = c("Focal adhesion", "ECM-receptor interaction", "Cell cycle", "DNA replication"),
    GeneRatio = c("5/50", "4/50", "6/50", "3/50"), BgRatio = c("200/8000", "90/8000", "150/8000", "36/8000"),
    pvalue = c(0.001, 0.002, 0.0005, 0.003), p.adjust = c(0.01, 0.02, 0.005, 0.03), qvalue = c(0.01, 0.02, 0.005, 0.03),
    # Focal adhesion and ECM-receptor share 4 of 5 genes (Jaccard 0.8); the other two are disjoint
    geneID = c("ITGA1/ITGB1/COL1A1/LAMB1/PTK2", "ITGA1/ITGB1/COL1A1/LAMB1", "CDK1/CCNB1/MCM2/CDC6/PCNA/MCM3", "MCM4/MCM5/POLA1"),
    Count = c(5L, 4L, 6L, 3L), stringsAsFactors = FALSE)
  rownames(df) <- df$ID
  new("enrichResult", result = df, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
      organism = "hsa", ontology = "KEGG", gene = character(), keytype = "SYMBOL",
      universe = character(), geneSets = list(), readable = FALSE)
}

test_that("non-GO results are reduced by Jaccard overlap, keeping the more significant term", {
  x <- suppressMessages(reduce_redundancy(make_kegg_enrich()))
  kept <- as.data.frame(x)$ID
  expect_setequal(kept, c("hsa04110", "hsa04510", "hsa03030"))
  expect_identical(attr(x, "reduction_method"), "jaccard")
  expect_identical(c(attr(x, "n_original"), attr(x, "n_reduced")), c(4L, 3L))
  # nothing is merged at a threshold above every overlap
  expect_identical(nrow(as.data.frame(suppressMessages(reduce_redundancy(make_kegg_enrich(), threshold = 0.9)))), 4L)
})

test_that("a compareClusterResult is reduced within each cluster, never across clusters", {
  e <- as.data.frame(make_kegg_enrich())
  df <- rbind(cbind(Cluster = "A", e), cbind(Cluster = "B", e[c(1, 3), ]))
  df$Cluster <- factor(df$Cluster)
  cc <- new("compareClusterResult", compareClusterResult = df, geneClusters = list(A = "x", B = "y"),
            fun = "enricher", .call = call("compareCluster"))
  red <- as.data.frame(suppressMessages(reduce_redundancy(cc)))
  expect_setequal(red$ID[red$Cluster == "A"], c("hsa04110", "hsa04510", "hsa03030"))
  expect_setequal(red$ID[red$Cluster == "B"], c("hsa04510", "hsa04110"))
})

test_that("GO results use Wang similarity, and the threshold controls the merge", {
  skip_if_not_installed("GO.db")
  df <- data.frame(ID = c("GO:0006915", "GO:0008283", "GO:0007049"),
                   Description = c("apoptotic process", "cell population proliferation", "cell cycle"),
                   GeneRatio = "2/5", BgRatio = "10/100", pvalue = c(0.001, 0.002, 0.003),
                   p.adjust = c(0.01, 0.02, 0.03), qvalue = 0.01, geneID = "A/B", Count = 2L, stringsAsFactors = FALSE)
  rownames(df) <- df$ID
  x <- new("enrichResult", result = df, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
           organism = "Homo sapiens", ontology = "BP", gene = character(), keytype = "SYMBOL",
           universe = character(), geneSets = list(), readable = FALSE)
  keep_all <- suppressMessages(reduce_redundancy(x, threshold = 1.01))
  keep_one <- suppressMessages(reduce_redundancy(x, threshold = -1))
  expect_identical(attr(keep_all, "reduction_method"), "Wang")
  expect_identical(nrow(as.data.frame(keep_all)), 3L)
  expect_identical(as.data.frame(keep_one)$ID, "GO:0006915")
})

test_that("n_pathways = NULL sends every significant term", {
  process <- getFromNamespace("process_enrichment_input", "clusterProfiler")
  x <- make_kegg_enrich()
  expect_identical(nrow(process(x, n_pathways = NULL)[[1]]$df), 4L)
  expect_identical(nrow(process(x, n_pathways = 2)[[1]]$df), 2L)
})

test_that("GSEA member genes (core_enrichment) reach the prompt", {
  fmt <- getFromNamespace(".format_pathway_text", "clusterProfiler")
  df <- data.frame(ID = "hsa04110", Description = "Cell cycle", NES = 2.1, p.adjust = 0.001,
                   core_enrichment = "CDK1/CCNB1/MCM2", stringsAsFactors = FALSE)
  expect_match(fmt(df), "CDK1/CCNB1/MCM2", fixed = TRUE)
})

test_that("the earlier argument names still switch the reduction", {
  rfp <- getFromNamespace(".reduce_for_prompt", "clusterProfiler")
  x <- make_kegg_enrich()
  expect_identical(nrow(as.data.frame(suppressMessages(rfp(x, reduce_go = FALSE)))), 4L)
  expect_identical(nrow(as.data.frame(suppressMessages(rfp(x, reduce = FALSE)))), 4L)
  expect_identical(nrow(as.data.frame(suppressMessages(rfp(x)))), 3L)
  # legacy calls keep their old default threshold of 0.7 (the 0.8 overlap still merges)
  expect_identical(attr(suppressMessages(rfp(x, reduce_go = TRUE)), "reduction_threshold"), 0.7)
})

test_that("all-NA padding rows neither reach the prompt nor change the similarity measure", {
  process <- getFromNamespace("process_enrichment_input", "clusterProfiler")
  x <- make_kegg_enrich()
  pad <- x@result[1, ]; pad[] <- NA; rownames(pad) <- "pad"
  x@result <- rbind(x@result, pad)
  expect_identical(nrow(process(x, n_pathways = NULL)[[1]]$df), 4L)
  red <- suppressMessages(reduce_redundancy(x))
  expect_identical(attr(red, "n_original"), 4L)
  expect_identical(attr(red, "reduction_method"), "jaccard")
})

test_that("an Entrez-keyed result with no symbol table reaches the prompt as gene symbols", {
  skip_if_not_installed("org.Hs.eg.db")
  rfp <- getFromNamespace(".readable_for_prompt", "clusterProfiler")
  df <- data.frame(ID = "GO:0007049", Description = "cell cycle", setSize = 3L, enrichmentScore = 0.8,
                   NES = 2.1, pvalue = 0.001, p.adjust = 0.01, qvalue = 0.01, rank = 1L,
                   leading_edge = "tags=100%", core_enrichment = "7157/4609/595", stringsAsFactors = FALSE)
  rownames(df) <- df$ID
  gl <- sort(c(`7157` = 3, `4609` = 2, `595` = 1, `1019` = -1), decreasing = TRUE)
  x <- new("gseaResult", result = df, organism = "Homo sapiens", setType = "BP", geneSets = list(),
           geneList = gl, keytype = "ENTREZID", permScores = matrix(), params = list(), readable = FALSE)
  y <- rfp(x)
  expect_true(y@readable)
  expect_setequal(strsplit(as.data.frame(y)$core_enrichment, "/")[[1]], c("TP53", "MYC", "CCND1"))
  # an object whose organism has no known OrgDb is left alone
  x@organism <- "Danio rerio"
  expect_false(rfp(x)@readable)
})
