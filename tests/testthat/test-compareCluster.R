context("Test compareCluster")

data(gcSample)
gc_small = gcSample[c("X1", "X2", 'X4')]
test_that("compareCluster + KEGG no formula", {
result = compareCluster(gc_small, fun="enrichKEGG",
                      organism="hsa", pvalueCutoff=0.05)
expect_is(result, 'compareClusterResult')
})


test_that("compareCluster + KEGG + formula", {
  gc_small_dat = do.call(rbind,
                         lapply(names(gc_small),
                                function(x) data.frame(X = x, entrezID = gc_small[[x]], stringsAsFactors = FALSE)
                                ))
  result = compareCluster(entrezID ~ X, data = gc_small_dat, fun="enrichKEGG",
                          organism="hsa", pvalueCutoff=0.05)
  expect_is(result, 'compareClusterResult')
})

test_that("compareCluster + GO no formula", {
  result = compareCluster(gc_small, fun="enrichGO", OrgDb = "org.Hs.eg.db",
                          pvalueCutoff=0.05)
  expect_is(result, 'compareClusterResult')
})

test_that("compareCluster + GO no formula, default value for function", {
  result = compareCluster(gc_small, OrgDb = "org.Hs.eg.db",
                          pvalueCutoff=0.05)
  expect_is(result, 'compareClusterResult')
})
