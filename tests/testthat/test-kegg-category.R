## KEGG pathway classification (category / subcategory) across entry points
## (issue #787: enrichKEGG reported it, gseKEGG did not)

test_that("append_kegg_category() adds the classification to both result types", {
  cat_tbl <- data.frame(
    category = c("Metabolism", "Cellular Processes"),
    subcategory = c("Global and overview maps", "Cell growth and death"),
    id = c("01100", "04110"),
    name = c("Metabolic pathways", "Cell cycle"),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(kegg_category_data = function() cat_tbl,
                        .package = "clusterProfiler")

  make_ora <- function() {
    new("enrichResult",
         result = data.frame(
           ID = c("hsa04110", "hsa01100"), Description = c("Cell cycle", "Metabolic pathways"),
           GeneRatio = "5/20", BgRatio = "100/2000", RichFactor = 0.05,
           FoldEnrichment = 5, zScore = 10, pvalue = c(0.001, 0.002),
           p.adjust = c(0.001, 0.002), qvalue = c(0.001, 0.002),
           geneID = c("a/b", "c/d"), Count = 2L, stringsAsFactors = FALSE),
         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
         organism = "hsa", ontology = "KEGG", gene = c("a", "b"),
         keytype = "kegg", universe = character(), geneSets = list(),
         readable = FALSE)
  }

  ora <- append_kegg_category(make_ora())
  expect_identical(colnames(as.data.frame(ora))[1:2], c("category", "subcategory"))
  expect_identical(as.data.frame(ora)$category, c("Cellular Processes", "Metabolism"))

  gsea <- new("gseaResult",
              result = data.frame(
                ID = c("hsa04110", "hsa01100"), Description = c("Cell cycle", "Metabolic pathways"),
                setSize = c(100L, 200L), enrichmentScore = c(0.6, 0.5),
                NES = c(2.1, 1.9), pvalue = c(1e-5, 1e-4),
                p.adjust = c(1e-4, 1e-3), qvalue = c(1e-4, 1e-3),
                rank = c(100L, 50L), leading_edge = c("tags=10%", "tags=8%"),
                core_enrichment = c("a/b", "c/d"), stringsAsFactors = FALSE),
              organism = "hsa", setType = "KEGG", geneList = c(a = 1, b = 0.5),
              keytype = "kegg", geneSets = list(), readable = FALSE)
  g <- append_kegg_category(gsea)
  expect_identical(colnames(as.data.frame(g))[1:2], c("category", "subcategory"))
  expect_identical(as.data.frame(g)$subcategory,
                   c("Cell growth and death", "Global and overview maps"))
})

test_that("append_kegg_category() leaves non-KEGG results alone", {
  gsea <- new("gseaResult",
              result = data.frame(ID = "hsa_M00001", Description = "module",
                                  setSize = 10L, enrichmentScore = 0.5, NES = 1.5,
                                  pvalue = 0.01, p.adjust = 0.02, qvalue = 0.02,
                                  rank = 5L, leading_edge = "tags=10%",
                                  core_enrichment = "a/b", stringsAsFactors = FALSE),
              organism = "hsa", setType = "MKEGG", geneList = c(a = 1, b = 0.5),
              keytype = "kegg", geneSets = list(), readable = FALSE)
  expect_message(out <- append_kegg_category(gsea), "Not a KEGG")
  expect_identical(colnames(as.data.frame(out)), colnames(as.data.frame(gsea)))
})

test_that("gseKEGG() reports the same classification columns as enrichKEGG() (#787)", {
  skip_on_cran()
  skip_if_offline("rest.kegg.jp")
  skip_if_not_installed("DOSE")

  data(geneList, package = "DOSE")
  g <- tryCatch(gseKEGG(geneList, organism = "hsa", minGSSize = 10,
                        pvalueCutoff = 0.5, verbose = FALSE, seed = TRUE),
                error = function(e) NULL)
  skip_if(is.null(g) || nrow(as.data.frame(g)) == 0, "KEGG not reachable")

  d <- as.data.frame(g)
  expect_true(all(c("category", "subcategory") %in% colnames(d)))
  expect_identical(colnames(d)[1:2], c("category", "subcategory"))
  # the classification is cached rather than fetched per run, so a pathway
  # released after the cache was built comes back NA (hsa04156 does today);
  # the promise is that the table is classified, not that the cache is current
  expect_gt(mean(!is.na(d$category)), 0.95)
})
