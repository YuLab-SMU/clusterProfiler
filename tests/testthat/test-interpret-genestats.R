library(testthat)

test_that("gene_stats tables from the common DE tools are normalised to gene / log2FC / padj", {
  gst <- getFromNamespace(".gene_stats_table", "clusterProfiler")
  limma <- data.frame(SYMBOL = c("A", "B"), logFC = c(1, -2), adj.P.Val = c(0.01, 0.2))
  deseq <- data.frame(row.names = c("A", "B"), log2FoldChange = c(1, -2), padj = c(0.01, 0.2))
  seurat <- data.frame(gene = c("A", "B"), avg_log2FC = c(1, -2), p_val_adj = c(0.01, 0.2))
  for (d in list(limma, deseq, seurat)) {
    st <- gst(d)
    expect_identical(names(st), c("gene", "log2FC", "padj"))
    expect_identical(st$gene, c("A", "B"))
    expect_equal(st$log2FC, c(1, -2))
  }
  # duplicates keep the most significant row; non-finite fold changes are dropped
  st <- gst(data.frame(gene = c("A", "B", "A", "C"), log2FC = c(1, -2, 3, NA), padj = c(0.01, 0.2, 0.001, 0.5)))
  expect_identical(st$gene, c("A", "B"))
  expect_equal(st$log2FC[st$gene == "A"], 3)
  expect_error(gst(data.frame(gene = "A", t = 1)), "log fold change column")
  expect_error(gst(list(gene = "A")), "data frame")
})

test_that("the gene-statistics section counts significant genes and lists each side", {
  gst <- getFromNamespace(".gene_stats_table", "clusterProfiler")
  gtx <- getFromNamespace(".gene_stats_text", "clusterProfiler")
  st <- gst(data.frame(gene = c("U1", "U2", "D1", "N1"), log2FC = c(2, 1, -3, 0.5), padj = c(0.001, 0.04, 0.002, 0.9)))
  txt <- gtx(st, n = 1)
  expect_match(txt, "4 genes measured; 2 up-regulated and 1 down-regulated at FDR < 0.05", fixed = TRUE)
  expect_match(txt, "up-regulated genes (up to 1):\nU1: log2FC 2.00, FDR 1.0e-03", fixed = TRUE)
  expect_match(txt, "down-regulated genes (up to 1):\nD1: log2FC -3.00", fixed = TRUE)
  expect_false(grepl("U2", txt))
  # without adjusted p-values: ranked by |log2FC|, no FDR claims
  txt0 <- gtx(gst(data.frame(row.names = c("X", "Y"), avg_log2FC = c(2, -1))))
  expect_false(grepl("FDR", txt0))
  expect_null(gtx(NULL))
})

test_that("term genes carry FDR only when per-gene statistics are supplied", {
  fct <- getFromNamespace(".get_fc_text", "clusterProfiler")
  df <- data.frame(ID = "T1", geneID = "A/B", stringsAsFactors = FALSE)
  fc <- c(A = 1.234, B = -2)
  expect_identical(fct(fc, df, NULL), "B:-2, A:1.23")
  expect_identical(fct(fc, df, NULL, padj = c(A = 0.001, B = NA)), "B:-2, A:1.23 (FDR 1.0e-03)")
})

test_that("the prompt keeps its old fold-change header when no label is given", {
  bds <- getFromNamespace(".build_data_sections", "clusterProfiler")
  p0 <- bds("T1", fold_change = "A:1")
  expect_match(p0, "Top Regulated Genes (Log2 Fold Change):\nA:1", fixed = TRUE)
  p1 <- bds("T1", fold_change = "A:1", fold_change_label = "Genes of the Enriched Terms (log2 fold change, FDR)",
            gene_stats = "2 genes measured.")
  expect_match(p1, "Genes of the Enriched Terms (log2 fold change, FDR):\nA:1", fixed = TRUE)
  expect_match(p1, "Gene-Level Statistics:\n2 genes measured.", fixed = TRUE)
})
