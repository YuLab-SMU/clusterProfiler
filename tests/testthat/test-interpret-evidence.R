library(testthat)

make_entrez_enrich <- function() {
  df <- data.frame(
    ID = c("GO:0006915", "GO:0008284"),
    Description = c("apoptotic process", "positive regulation of proliferation"),
    GeneRatio = c("2/5", "3/5"), BgRatio = c("10/100", "20/100"),
    pvalue = c(0.001, 0.002), p.adjust = c(0.01, 0.02), qvalue = c(0.01, 0.02),
    geneID = c("TP53/BAX", "MYC/CCND1/CDK4"), Count = c(2L, 3L),
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$ID
  new("enrichResult", result = df, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
      organism = "Homo sapiens", ontology = "BP", gene = c("7157", "581", "4609", "595", "1019"),
      keytype = "ENTREZID", universe = character(), geneSets = list(), readable = TRUE,
      gene2Symbol = c(`7157` = "TP53", `581` = "BAX", `4609` = "MYC", `595` = "CCND1", `1019` = "CDK4"))
}

test_that("input gene IDs reach the prompt as symbols, and unmapped IDs pass through", {
  ids_to_symbols <- getFromNamespace(".ids_to_symbols", "clusterProfiler")
  x <- make_entrez_enrich()
  expect_identical(ids_to_symbols(c("7157", "4609"), x), c("TP53", "MYC"))
  expect_identical(ids_to_symbols(c("7157", "999999"), x), c("TP53", "999999"))
  # an object with no symbol table is left unchanged
  x@gene2Symbol <- character()
  expect_identical(ids_to_symbols(c("7157", "4609"), x), c("7157", "4609"))
})

test_that("gene_fold_change is renamed with the same table so it still matches", {
  fc_to_symbols <- getFromNamespace(".fc_to_symbols", "clusterProfiler")
  x <- make_entrez_enrich()
  fc <- c(`7157` = 2.5, `4609` = 1.2, `999999` = 0.4)
  expect_identical(fc_to_symbols(fc, x), c(TP53 = 2.5, MYC = 1.2, `999999` = 0.4))
  expect_null(fc_to_symbols(NULL, x))
})

test_that("the cell-type prompt lists marker symbols, not Entrez numbers", {
  process <- getFromNamespace("process_enrichment_input", "clusterProfiler")
  res <- process(make_entrez_enrich(), n_pathways = 20)
  genes <- res[[1]]$genes
  expect_true(all(c("TP53", "MYC") %in% genes))
  expect_false(any(grepl("^[0-9]+$", genes)))
})

test_that("evidence_status labels every named gene as observed or out of evidence", {
  x <- make_entrez_enrich()
  rec <- structure(list(
    overview = "TP53 and BAX drive apoptosis; MYC is induced. AKT1 may act upstream.",
    hypothesis = "Knockdown of MTOR would reduce CDK4."
  ), class = c("interpretation", "list"))
  st <- evidence_status(rec, x, symbols = c("TP53", "BAX", "MYC", "CDK4", "AKT1", "MTOR"))
  expect_setequal(st$gene[st$status == "observed"], c("TP53", "BAX", "MYC", "CDK4"))
  expect_setequal(st$gene[st$status == "out_of_evidence"], c("AKT1", "MTOR"))
  expect_setequal(names(st), c("cluster", "field", "gene", "status"))
})

test_that("evidence_status ignores capitalised words unless they are symbols or evidence", {
  x <- make_entrez_enrich()
  rec <- structure(list(overview = "DNA damage and ATP levels, with TP53. GO BP terms. Cell cycle."),
                   class = c("interpretation", "list"))
  st <- evidence_status(rec, x, symbols = c("TP53", "ATP", "DNA"))
  expect_identical(st$gene, "TP53")
  # without a vocabulary only evidence genes are recognised
  st2 <- evidence_status(rec, x)
  expect_identical(st2$gene, "TP53")
})

test_that("evidence_status reads compareCluster records against their own cluster", {
  df <- data.frame(Cluster = c("A", "B"), ID = c("GO:1", "GO:2"), Description = c("a", "b"),
                   geneID = c("TP53/BAX", "MYC"), p.adjust = c(0.01, 0.01), stringsAsFactors = FALSE)
  recs <- structure(list(A = list(overview = "TP53 and MYC"), B = list(overview = "MYC and TP53")),
                    class = c("interpretation_list", "list"))
  st <- evidence_status(recs, df, symbols = c("TP53", "BAX", "MYC"))
  expect_identical(st$status[st$cluster == "A" & st$gene == "MYC"], "out_of_evidence")
  expect_identical(st$status[st$cluster == "B" & st$gene == "MYC"], "observed")
  expect_identical(st$status[st$cluster == "B" & st$gene == "TP53"], "out_of_evidence")
})

test_that("evidence_status returns an empty frame when no gene is named", {
  rec <- structure(list(overview = "Nothing specific."), class = c("interpretation", "list"))
  st <- evidence_status(rec, make_entrez_enrich())
  expect_equal(nrow(st), 0)
  expect_setequal(names(st), c("cluster", "field", "gene", "status"))
})
