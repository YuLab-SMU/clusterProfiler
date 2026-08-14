library(clusterProfiler)
library(testthat)

test_that("GSEA wrappers expose eps with the enrichit default", {
    expect_equal(formals(clusterProfiler::GSEA)$eps, 1e-10)
    expect_equal(formals(clusterProfiler::gseGO)$eps, 1e-10)
    expect_equal(formals(clusterProfiler::gseMKEGG)$eps, 1e-10)
    expect_equal(formals(clusterProfiler::gseKEGG)$eps, 1e-10)
})

test_that("GSEA wrappers expose a seed argument defaulting to FALSE", {
    expect_equal(formals(clusterProfiler::GSEA)$seed, FALSE)
    expect_equal(formals(clusterProfiler::gseGO)$seed, FALSE)
    expect_equal(formals(clusterProfiler::gseMKEGG)$seed, FALSE)
    expect_equal(formals(clusterProfiler::gseKEGG)$seed, FALSE)
})

test_that("GSEA forwards eps and extra enrichit arguments", {
    captured <- NULL

    local_mocked_bindings(
        gsea_gson = function(...) {
            captured <<- list(...)
            NULL
        },
        .package = "enrichit"
    )

    geneList <- c(geneA = 2, geneB = 1, geneC = -1)
    term2gene <- data.frame(
        term = c("set1", "set1", "set2"),
        gene = c("geneA", "geneB", "geneC")
    )

    clusterProfiler::GSEA(
        geneList = geneList,
        TERM2GENE = term2gene,
        eps = 0,
        sampleSize = 201,
        seed = 42,
        verbose = FALSE
    )

    expect_equal(captured$eps, 0)
    expect_equal(captured$sampleSize, 201)
    expect_equal(captured$seed, 42)
    expect_s4_class(captured$gson, "GSON")
})
