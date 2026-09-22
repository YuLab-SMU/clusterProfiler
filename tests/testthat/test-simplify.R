## simplify() on a GSEA result whose gene sets are GO terms (#753)

test_that("infer_go_ontology() derives the ontology from GO IDs", {
  skip_if_not_installed("GO.db")

  # both are biological process terms
  expect_equal(infer_go_ontology(c("GO:0008150", "GO:0009987")), "BP")

  # a BP term plus a cellular component term spans ontologies
  expect_equal(infer_go_ontology(c("GO:0008150", "GO:0005575")), "GOALL")

  # not GO terms at all -> empty string, so the caller keeps its own error
  expect_equal(infer_go_ontology(c("path:hsa00010", "path:hsa00020")), "")

  # mixed GO and non-GO is not a GO collection either
  expect_equal(infer_go_ontology(c("GO:0008150", "hsa00010")), "")

  # empty / NA input must not error
  expect_equal(infer_go_ontology(character(0)), "")
  expect_equal(infer_go_ontology(NA_character_), "")
})

test_that("simplify() works on a GSEA result over a GO collection (#753)", {
  skip_if_not_installed("GO.db")
  skip_if_not_installed("org.Hs.eg.db")

  # a small GO BP collection built from real annotations; the OrgDb object is
  # fetched from the suggested package without attaching it
  odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  keys_all <- head(AnnotationDbi::keys(odb, "ENTREZID"), 800)
  go2gene <- suppressMessages(
    AnnotationDbi::select(odb, keys = keys_all,
                          columns = "GO", keytype = "ENTREZID")
  )
  go2gene <- stats::na.omit(go2gene)
  ont <- AnnotationDbi::Ontology(GO.db::GOTERM)[go2gene$GO]
  bp <- go2gene[!is.na(ont) & ont == "BP", c("GO", "ENTREZID")]
  skip_if(nrow(bp) < 20, "not enough GO annotations available")

  geneList <- sort(setNames(rnorm(length(unique(bp$ENTREZID))), unique(bp$ENTREZID)),
                   decreasing = TRUE)

  x <- GSEA(geneList, TERM2GENE = bp, pvalueCutoff = 1, minGSSize = 2,
            maxGSSize = 1000, verbose = FALSE)
  skip_if(is.null(x), "GSEA returned no result")

  # @setType carries the collection name, not an ontology -> simplify() must
  # still work by inferring the ontology from the GO IDs
  expect_false(x@setType %in% c("BP", "MF", "CC", "GOALL"))

  y <- simplify(x, cutoff = 0.7, measure = "Wang")
  expect_s4_class(y, "gseaResult")
  expect_true(nrow(as.data.frame(y)) <= nrow(as.data.frame(x)))
})

test_that("simplify() still refuses non-GO gene sets with a clear error (#753)", {
  geneList <- sort(setNames(c(2, 1, -1, -2), c("g1", "g2", "g3", "g4")), decreasing = TRUE)
  t2g <- data.frame(
    term = c("path:hsa00010", "path:hsa00010", "path:hsa00020", "path:hsa00020"),
    gene = c("g1", "g2", "g3", "g4")
  )
  x <- GSEA(geneList, TERM2GENE = t2g, pvalueCutoff = 1, minGSSize = 1, verbose = FALSE)
  skip_if(is.null(x), "GSEA returned no result")

  expect_error(simplify(x, cutoff = 0.7), "gene sets are GO terms")
})
