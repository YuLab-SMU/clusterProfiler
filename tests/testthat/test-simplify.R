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

test_that("simplify() works on an enricher() result over a GO collection (#369)", {
  skip_if_not_installed("GO.db")
  skip_if_not_installed("org.Hs.eg.db")

  odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  keys_all <- head(AnnotationDbi::keys(odb, "ENTREZID"), 1200)
  go2gene <- suppressMessages(
    AnnotationDbi::select(odb, keys = keys_all, columns = "GO", keytype = "ENTREZID")
  )
  go2gene <- stats::na.omit(go2gene)
  ont <- AnnotationDbi::Ontology(GO.db::GOTERM)[go2gene$GO]
  bp <- go2gene[!is.na(ont) & ont == "BP", c("GO", "ENTREZID")]
  skip_if(nrow(bp) < 20, "not enough GO annotations available")

  # a subset of the annotated genes, so that some terms are enriched at all
  set.seed(1)
  gene <- sample(unique(bp$ENTREZID), min(150, length(unique(bp$ENTREZID))))
  res <- enricher(gene, TERM2GENE = bp, pvalueCutoff = 1, minGSSize = 5)
  skip_if(is.null(res) || nrow(as.data.frame(res)) == 0, "enricher returned no result")

  # @ontology is not a GO ontology here, yet simplify() must work without the
  # user assigning the slot by hand
  expect_false(res@ontology %in% c("BP", "MF", "CC", "GOALL"))

  y <- simplify(res, cutoff = 0.7, measure = "Wang")
  expect_s4_class(y, "enrichResult")
  expect_true(nrow(as.data.frame(y)) <= nrow(as.data.frame(res)))
})

test_that("simplify() still refuses a non-GO enricher() result (#369)", {
  # gene list covers only part of each set, so terms are tested and returned
  t2g <- data.frame(term = rep(paste0("set", 1:5), each = 40),
                    gene = paste0("g", 1:200))
  gene <- paste0("g", 1:60)
  res <- enricher(gene, TERM2GENE = t2g, pvalueCutoff = 1, minGSSize = 2,
                  universe = paste0("g", 1:400))
  skip_if(is.null(res) || nrow(as.data.frame(res)) == 0, "enricher returned no result")

  expect_error(simplify(res, cutoff = 0.7), "gene sets are GO terms")
})

## Gene-overlap redundancy removal for non-GO collections (#743)

make_overlap_result <- function() {
  # three overlapping sets plus one unrelated set, all tested on 200 genes
  t2g <- data.frame(
    term = c(rep("A", 40), rep("B", 40), rep("C", 40), rep("D", 40)),
    gene = c(paste0("g", 1:40),             # A
             paste0("g", 21:60),            # B overlaps A by 20 / union 60
             paste0("g", 1:40),             # C is identical to A
             paste0("u", 1:40))             # D is unrelated
  )
  gene <- c(paste0("g", 1:30), "u1", "u2")
  enricher(gene, TERM2GENE = t2g, pvalueCutoff = 1, qvalueCutoff = 1,
           minGSSize = 2, universe = c(paste0("g", 1:60), paste0("u", 1:40)))
}

test_that("simplify(measure = 'JC') removes redundant gene sets (#743)", {
  res <- make_overlap_result()
  skip_if(is.null(res) || nrow(as.data.frame(res)) == 0, "enricher returned no result")

  y <- simplify(res, cutoff = 0.7, measure = "JC")
  keep <- as.data.frame(y)$ID

  # A and C are the same gene set, so only one of them survives
  expect_equal(sum(c("A", "C") %in% keep), 1)
  # B shares 20 of 60 genes with A: Jaccard 1/3 < 0.7, so B is kept
  expect_true("B" %in% keep)
  # D does not overlap anything
  expect_true("D" %in% keep)
  expect_s4_class(y, "enrichResult")

  # A vs B is only 20/60 = 0.33, so a 0.5 cutoff still keeps B while
  # collapsing the identical pair A / C
  y2 <- simplify(res, cutoff = 0.5, measure = "JC")
  keep2 <- as.data.frame(y2)$ID
  expect_equal(sum(c("A", "C") %in% keep2), 1)
  expect_true("B" %in% keep2)
})

test_that("simplify(measure = 'overlap') uses shared/smaller-set (#743)", {
  res <- make_overlap_result()
  skip_if(is.null(res) || nrow(as.data.frame(res)) == 0, "enricher returned no result")

  # A vs B: 20 shared of 40 each -> 0.5, above a 0.4 cutoff
  y <- simplify(res, cutoff = 0.4, measure = "overlap")
  keep <- as.data.frame(y)$ID
  expect_equal(sum(c("A", "B", "C") %in% keep), 1)
})

test_that("simplify(measure = 'JC') keeps the representative set by `by`", {
  res <- make_overlap_result()
  skip_if(is.null(res) || nrow(as.data.frame(res)) == 0, "enricher returned no result")

  # A and C are identical gene sets but differ in p.adjust, so the better one
  # must be the survivor no matter which row it sits in
  df <- as.data.frame(res)
  df$p.adjust[df$ID == "A"] <- 0.001
  df$p.adjust[df$ID == "C"] <- 0.5
  res@result <- df

  keep <- as.data.frame(simplify(res, cutoff = 0.7, measure = "JC"))$ID
  expect_true("A" %in% keep)
  expect_false("C" %in% keep)
})

test_that("simplify(measure = 'JC') works on a KEGG result (#743)", {
  skip_on_cran()
  skip_if_offline("rest.kegg.jp")

  data(geneList, package = "DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  kk <- tryCatch(
    enrichKEGG(gene, organism = "hsa", pvalueCutoff = 0.5),
    error = function(e) NULL
  )
  skip_if(is.null(kk) || nrow(as.data.frame(kk)) < 2, "KEGG not reachable")

  y <- simplify(kk, cutoff = 0.7, measure = "JC")
  expect_s4_class(y, "enrichResult")
  expect_true(nrow(as.data.frame(y)) <= nrow(as.data.frame(kk)))
  expect_true(nrow(as.data.frame(y)) > 0)
})

test_that("simplify() points non-GO users to measure = 'JC' (#743)", {
  t2g <- data.frame(term = rep(paste0("set", 1:3), each = 40),
                    gene = paste0("g", 1:120))
  gene <- paste0("g", 1:50)
  res <- enricher(gene, TERM2GENE = t2g, pvalueCutoff = 1, minGSSize = 2,
                  universe = paste0("g", 1:300))
  skip_if(is.null(res) || nrow(as.data.frame(res)) == 0, "enricher returned no result")

  expect_error(simplify(res, cutoff = 0.7), 'measure = "JC"')
})

test_that("gene-overlap simplify does not delete a term whose only neighbour is gone (#743)", {
  # A~B = 0.739 and B~C = 0.739 are both above 0.7, but A~C = 0.538 is not.
  # A absorbs B; C's only similar partner was B, which is already removed, so
  # C must survive. The old grouping tested "is anything similar to C" and
  # dropped C behind the removed B.
  A <- paste0("g", 1:40)
  B <- c(paste0("g", 1:34), paste0("g", 100:105))
  C <- c(paste0("g", 1:28), paste0("g", 100:105), paste0("g", 200:205))
  t2g <- data.frame(term = rep(c("A", "B", "C"), each = 40), gene = c(A, B, C))
  res <- enricher(paste0("g", 1:20), TERM2GENE = t2g, pvalueCutoff = 1,
                  qvalueCutoff = 1, minGSSize = 2, universe = paste0("g", 1:400))
  skip_if(is.null(res) || nrow(as.data.frame(res)) < 3, "enricher returned too few terms")

  df <- as.data.frame(res)
  df$p.adjust <- c(0.01, 0.02, 0.03)[match(df$ID, c("A", "B", "C"))]
  res@result <- df

  sim <- clusterProfiler:::gene_overlap_similarity(
    clusterProfiler:::gene_overlap_sets(res, df), "JC")
  skip_if(unname(sim[1, 2]) <= 0.7 || unname(sim[2, 3]) <= 0.7 ||
            unname(sim[1, 3]) > 0.7, "fixture no longer has the transitivity shape")

  keep <- as.data.frame(simplify(res, cutoff = 0.7, measure = "JC"))$ID
  expect_true("A" %in% keep)   # best p.adjust
  expect_false("B" %in% keep)  # absorbed by A
  expect_true("C" %in% keep)   # only neighbour B is gone, so it must stay
})

test_that("gene-overlap simplify honours select_fun (#743)", {
  # A and B are the same gene set, so exactly one survives, and which one is
  # decided by `select_fun` -- the implementation used to sort ascending
  # regardless, so select_fun = max still returned the min-p term
  t2g <- data.frame(term = rep(c("A", "B"), each = 10), gene = rep(paste0("g", 1:10), 2))
  res <- enricher(paste0("g", 1:5), TERM2GENE = t2g, pvalueCutoff = 1,
                  qvalueCutoff = 1, minGSSize = 2, universe = paste0("g", 1:200))
  skip_if(is.null(res) || nrow(as.data.frame(res)) < 2, "enricher returned too few terms")

  df <- as.data.frame(res)
  df$p.adjust <- c(0.01, 0.9)[match(df$ID, c("A", "B"))]
  res@result <- df

  expect_identical(as.data.frame(simplify(res, 0.7, measure = "JC"))$ID, "A")
  expect_identical(
    as.data.frame(simplify(res, 0.7, measure = "JC", select_fun = max))$ID, "B")
})

test_that("gene-overlap simplify breaks ties towards the larger gene set (#743)", {
  # BIG and SML share 10 of BIG's 12 genes (Jaccard 0.83 > 0.7) and have equal
  # p.adjust, so the larger set must be the representative even though SML is
  # listed first
  t2g <- data.frame(term = c(rep("BIG", 12), rep("SML", 10)),
                    gene = c(paste0("g", 1:12), paste0("g", 1:10)))
  res <- enricher(paste0("g", 1:5), TERM2GENE = t2g, pvalueCutoff = 1,
                  qvalueCutoff = 1, minGSSize = 2, universe = paste0("g", 1:200))
  skip_if(is.null(res) || nrow(as.data.frame(res)) < 2, "enricher returned too few terms")

  df <- as.data.frame(res)
  df$p.adjust <- 0.01
  res@result <- df
  skip_if(!setequal(df$ID, c("BIG", "SML")), "fixture no longer has both sets")

  keep <- as.data.frame(simplify(res, 0.7, measure = "JC"))$ID
  expect_identical(keep, "BIG")
})

test_that("compareCluster simplify(measure = 'JC') reduces within clusters only", {
  genes <- paste0("g", 1:10)
  t2g <- data.frame(term = rep(c("A", "B", "C"), each = 10), gene = rep(genes, 3))
  cc <- compareCluster(list(c1 = paste0("g", 1:5), c2 = paste0("g", 1:5)),
                       fun = "enricher", TERM2GENE = t2g, pvalueCutoff = 1,
                       qvalueCutoff = 1, minGSSize = 2, universe = paste0("g", 1:200))
  skip_if(is.null(cc) || nrow(as.data.frame(cc)) == 0, "compareCluster returned nothing")

  # identical clusters: each keeps one representative of the A/B pair
  s <- simplify(cc, 0.7, measure = "JC")
  d <- as.data.frame(s)
  expect_equal(nrow(d), 2L)
  expect_setequal(d$Cluster, c("c1", "c2"))
  expect_true(all(d$ID == "A"))  # best p.adjust per cluster

  # clusters are independent: C is identical to A but sits only in c2, and
  # c1's A must not remove it
  raw <- as.data.frame(cc)
  keep_rows <- (raw$Cluster == "c1" & raw$ID %in% c("A", "B")) |
    (raw$Cluster == "c2" & raw$ID == "C")
  cc2 <- cc
  cc2@compareClusterResult <- raw[keep_rows, ]
  s2 <- simplify(cc2, 0.7, measure = "JC")
  expect_setequal(as.data.frame(s2)$ID, c("A", "C"))
})
