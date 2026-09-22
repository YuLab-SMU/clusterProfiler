# data(geneList, package="DOSE")

# mydf <- data.frame(Entrez=names(geneList), FC=geneList)
# mydf <- mydf[abs(mydf$FC) > 1,]
# mydf$group <- "upregulated"
# mydf$group[mydf$FC < 0] <- "downregulated"
# mydf$othergroup <- "A"
# mydf$othergroup[abs(mydf$FC) > 2] <- "B"

# test_that("enrichGO formula interface works", {
#   formula_res <- compareCluster(
#     Entrez ~ group + othergroup,
#     data = mydf, fun = "enrichGO", OrgDb = org.Hs.eg.db)

#   expect_true(is(formula_res, "compareClusterResult"))
#   expect_equal(formula_res@fun, "enrichGO")
#   expect_true(all(sapply(formula_res@geneClusters, function(x) is.character(x))))
# })

# test_that("gseGO formula interface works", {
#   formula_res <- compareCluster(
#     Entrez | FC ~ group + othergroup,
#     data = mydf, fun = "gseGO", OrgDb = org.Hs.eg.db)

#   expect_true(is(formula_res, "compareClusterResult"))
#   expect_equal(formula_res@fun, "gseGO")
#   expect_true(all(sapply(formula_res@geneClusters, function(x) is.numeric(x))))
# })
test_that("a non-character `universe` warns instead of being silently dropped (#654)", {
  # compareCluster() wraps each per-cluster call in suppressMessages(), which used
  # to hide enrichGO()'s "universe ... will be ignored" notice completely.
  skip_if_not_installed("org.Hs.eg.db")
  gene <- c("1", "2", "3")
  expect_warning(
    compareCluster(
      list(A = gene, B = gene),
      fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP",
      universe = as.numeric(c(gene, "4", "5"))
    ),
    "universe"
  )
})

test_that("a character `universe` does not raise the universe warning (#654)", {
  skip_if_not_installed("org.Hs.eg.db")
  gene <- c("1", "2", "3")
  seen <- character()
  withCallingHandlers(
    try(
      compareCluster(
        list(A = gene, B = gene),
        fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP",
        universe = c(gene, "4", "5")
      ),
      silent = TRUE
    ),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # "No enrichment found ..." is expected for these dummy IDs; the universe
  # warning is what must NOT appear.
  expect_false(any(grepl("universe", seen, ignore.case = TRUE)))
})
