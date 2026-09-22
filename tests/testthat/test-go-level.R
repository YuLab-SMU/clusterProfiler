## GO level information: go_level_map() / getGOLevel() / add_go_level() (#793)

small_go_result <- function(ont = "BP", n = 200) {
    odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    gene <- head(AnnotationDbi::keys(odb, "ENTREZID"), n)
    enrichGO(gene, OrgDb = "org.Hs.eg.db", ont = ont,
             pvalueCutoff = 1, qvalueCutoff = 1)
}

test_that("go_level_map() places the ontology root at level 1", {
    skip_if_not_installed("GO.db")

    for (ont in c("BP", "CC", "MF")) {
        m <- clusterProfiler:::go_level_map(ont)
        root <- c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")[[ont]]

        expect_equal(unname(m[[root]]), 1L)
        expect_true(all(m >= 1L))
        expect_true(length(m) > 100)
    }

    expect_error(clusterProfiler:::go_level_map("nope"), "ontology")
})

test_that("getGOLevel() returns the frontier at the requested level", {
    skip_if_not_installed("GO.db")

    expect_equal(clusterProfiler:::getGOLevel("BP", 1), "GO:0008150")
    expect_equal(clusterProfiler:::getGOLevel("CC", 1), "GO:0005575")

    # asking for several levels gives their union
    u <- clusterProfiler:::getGOLevel("BP", c(2, 3))
    expect_equal(
        length(u),
        length(clusterProfiler:::getGOLevel("BP", 2)) +
            length(clusterProfiler:::getGOLevel("BP", 3))
    )

    # a level beyond the DAG depth has no terms
    expect_equal(length(clusterProfiler:::getGOLevel("BP", 40)), 0)
})

test_that("add_go_level() labels every term with its ontology level (#793)", {
    skip_if_not_installed("GO.db")
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")

    x <- small_go_result("BP")
    skip_if(nrow(as.data.frame(x)) == 0, "no enriched terms")

    y <- add_go_level(x)
    d <- as.data.frame(y)

    expect_s4_class(y, "enrichResult")
    expect_true("level" %in% names(d))
    expect_false(any(is.na(d$level)))
    expect_true(all(d$level >= 1))

    # each reported level agrees with the frontier definition it came from
    for (i in seq_len(min(5, nrow(d)))) {
        ont <- as.character(AnnotationDbi::Ontology(GO.db::GOTERM)[d$ID[i]])
        expect_true(d$ID[i] %in% clusterProfiler:::getGOLevel(ont, d$level[i]))
    }

    # filtering by level is what the column is for
    band <- subset(d, level >= 3 & level <= 6)
    expect_true(all(band$level >= 3 & band$level <= 6))
})

test_that("add_go_level() handles a multi-ontology result term by term", {
    skip_if_not_installed("GO.db")
    skip_if_not_installed("org.Hs.eg.db")

    x <- small_go_result("ALL")
    skip_if(nrow(as.data.frame(x)) == 0, "no enriched terms")

    d <- as.data.frame(add_go_level(x))
    expect_false(any(is.na(d$level)))

    # a CC term must be levelled in CC, not in BP
    ont_of <- as.character(AnnotationDbi::Ontology(GO.db::GOTERM)[d$ID])
    expect_true(length(unique(ont_of)) > 1)
    for (i in seq_len(min(5, nrow(d)))) {
        expect_true(d$ID[i] %in% clusterProfiler:::getGOLevel(ont_of[i], d$level[i]))
    }
})

test_that("gofilter() and dropGO() still select the requested levels", {
    skip_if_not_installed("GO.db")
    skip_if_not_installed("org.Hs.eg.db")

    x <- small_go_result("BP")
    skip_if(nrow(as.data.frame(x)) == 0, "no enriched terms")

    filtered <- gofilter(x, level = 3)
    expect_true(all(as.data.frame(filtered)$ID %in%
        clusterProfiler:::getGOLevel("BP", 3)))

    dropped <- dropGO(x, level = 3)
    expect_false(any(as.data.frame(dropped)$ID %in%
        clusterProfiler:::getGOLevel("BP", 3)))

    # dropping and keeping the same level are complementary
    expect_equal(
        nrow(as.data.frame(filtered)) + nrow(as.data.frame(dropped)),
        nrow(as.data.frame(x))
    )
})
