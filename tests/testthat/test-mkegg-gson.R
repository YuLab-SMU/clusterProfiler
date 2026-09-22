## KEGG Module analysis from a GSON, i.e. without contacting KEGG (#623)

mkegg_gson <- function() {
    gene_sets <- list(
        M00001 = c("g1", "g2", "g3"),
        M00002 = c("g2", "g3", "g4"),
        M00003 = c("g5", "g6")
    )
    gson::gson(
        gsid2gene = data.frame(
            gsid = rep(names(gene_sets), lengths(gene_sets)),
            gene = unlist(gene_sets)
        ),
        gsid2name = data.frame(
            gsid = names(gene_sets),
            name = paste("module", names(gene_sets))
        ),
        species = "test", gsname = "MKEGG", keytype = "kegg",
        version = "test", accessed_date = as.character(Sys.Date())
    )
}

test_that("enrichMKEGG() accepts a GSON instead of a species name (#623)", {
    skip_if_not_installed("gson")
    x <- enrichMKEGG(c("g1", "g2", "g3", "g4"), organism = mkegg_gson(),
                     pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)

    expect_s4_class(x, "enrichResult")
    expect_true(nrow(as.data.frame(x)) > 0)
    # the GSON carries the metadata a species code cannot
    expect_equal(x@organism, "test")
    expect_equal(x@ontology, "MKEGG")
    expect_equal(x@keytype, "kegg")
})

test_that("gseMKEGG() accepts a GSON instead of a species name (#623)", {
    skip_if_not_installed("gson")
    geneList <- sort(setNames(c(3, 2, 1, -1, -2, -3),
                              c("g1", "g2", "g3", "g4", "g5", "g6")),
                     decreasing = TRUE)
    x <- gseMKEGG(geneList = geneList, organism = mkegg_gson(),
                  pvalueCutoff = 1, minGSSize = 1, verbose = FALSE)

    expect_s4_class(x, "gseaResult")
    expect_equal(x@organism, "test")
    expect_equal(x@setType, "MKEGG")
    expect_equal(x@keytype, "kegg")
})

test_that("enrichMKEGG()/gseMKEGG() reject an organism that is neither", {
    expect_error(enrichMKEGG(c("g1"), organism = 42),
                 "species name or a GSON object")
    expect_error(gseMKEGG(geneList = c(a = 1), organism = list()),
                 "species name or a GSON object")
})
