data(geneList, package="DOSE")

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

test_that("enrichGO formula interface works", {
  formula_res <- compareCluster(
    Entrez ~ group + othergroup,
    data = mydf, fun = "enrichGO", OrgDb = org.Hs.eg.db)

  expect_true(is(formula_res, "compareClusterResult"))
  expect_equal(formula_res@fun, "enrichGO")
  expect_true(all(sapply(formula_res@geneClusters, function(x) is.character(x))))
})

test_that("gseGO formula interface works", {
  formula_res <- compareCluster(
    Entrez | FC ~ group + othergroup,
    data = mydf, fun = "gseGO", OrgDb = org.Hs.eg.db)

  expect_true(is(formula_res, "compareClusterResult"))
  expect_equal(formula_res@fun, "gseGO")
  expect_true(all(sapply(formula_res@geneClusters, function(x) is.numeric(x))))
})