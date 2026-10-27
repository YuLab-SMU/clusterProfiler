library(clusterProfiler)

context("enrichGO")

test_that("non-ENTREZID universe is converted together with gene", {
    skip_if_not_installed("org.Hs.eg.db")

    raw_accnum_gene <- c(
        "O00411", "Q9BQP7", "Q8TCS8", "Q8WV74", "O14734", "Q96C36",
        "P48728", "Q5T440", "Q9NYK5", "Q92665", "Q92947", "Q969S9",
        "Q8WWH5", "Q9BSH4", "M0QWZ7", "P09543", "A0A0G2JL52", "Q7L8L6",
        "Q92506", "Q96GW9", "A0A0C4DGA2", "Q9BTZ2", "Q9NUB1", "P35914",
        "F6T1Q0", "P14174", "Q15067", "P30085", "J3KS15", "Q9BY49",
        "Q96GK7", "Q08426", "P52948", "A0A182DWF2", "O00161", "Q9P0Z9",
        "P22570", "P62760", "P43155", "Q99714", "P43897", "O75439",
        "Q96I99", "Q86X76", "F5H5I6", "Q8NE62", "H7BXY3", "Q9NVV4",
        "A0A0B4J1R2", "Q02252", "Q9BQ69", "P23378", "Q5JTZ9", "C9JQS9",
        "P04179", "P51649", "G3V325", "A0A2R8Y6Y7", "C9JRZ6", "P51659",
        "P36551", "Q9HCC0", "P04040", "P42765", "O15118", "F5GX62",
        "Q9NPJ3", "P55084", "O14561", "Q13011", "P35270", "P36957",
        "F6TLX2", "R4GN18", "Q9HBH1", "Q9Y624", "P53999", "Q9UIJ7",
        "Q8N5N7", "P24752", "Q6YN16", "P30084", "P23141", "P11413",
        "Q9BQA1", "H0YFD6", "P42126", "Q9UMS0", "Q13268", "Q00059",
        "Q86TX2", "P49748", "Q9HAV7", "O00330", "P61604", "P00390",
        "Q04837", "O76021", "P00352", "P27695", "Q9Y4W6", "P34897",
        "P28838", "P45954", "Q8TD30"
    )

    accnum_map <- suppressWarnings(bitr(
        raw_accnum_gene,
        fromType = "ACCNUM",
        toType = "ENTREZID",
        OrgDb = "org.Hs.eg.db"
    ))
    accnum_gene <- unique(accnum_map$ACCNUM)
    entrez_gene <- unique(accnum_map$ENTREZID)

    accnum_res <- enrichGO(
        accnum_gene,
        OrgDb = "org.Hs.eg.db",
        keyType = "ACCNUM",
        ont = "BP",
        universe = accnum_gene,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 1
    )
    expect_s4_class(accnum_res, "enrichResult")
    expect_identical(accnum_res@keytype, "ACCNUM")

    entrez_res <- enrichGO(
        entrez_gene,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENTREZID",
        ont = "BP",
        universe = entrez_gene,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 1
    )
    expect_s4_class(entrez_res, "enrichResult")

    compare_cols <- c(
        "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
        "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "Count"
    )

    expect_equal(
        as.data.frame(accnum_res)[, compare_cols],
        as.data.frame(entrez_res)[, compare_cols],
        tolerance = 1e-12
    )
})

test_that("the GO cache is not reused for an unidentifiable OrgDb (#232)", {
    ## The cache holds the all-ontology annotation and `ont` is not part of
    ## the key, so a cache written by any ont serves any other ont. Keying on
    ## `ont` used to evict the annotation on every ontology switch:
    ## five alternating BP/MF enrichGO() calls took 85.5 s and 11.9 s after
    ## the key change (human SYMBOL, 2026-10-27).
    mk_env <- function(organism = "Homo sapiens", keytype = "ENTREZID",
                       goAnno = TRUE) {
        e <- new.env()
        assign("organism", organism, envir = e)
        assign("keytype", keytype, envir = e)
        if (goAnno) assign("goAnno", data.frame(a = 1), envir = e)
        e
    }

    # a matching request can reuse the cache
    expect_true(go_cache_usable(mk_env(), "Homo sapiens", "BP", "ENTREZID"))
    # a cache written under one ont serves a request for any other ont
    expect_true(go_cache_usable(mk_env(), "Homo sapiens", "MF", "ENTREZID"))
    expect_true(go_cache_usable(mk_env(), "Homo sapiens", "ALL", "ENTREZID"))
    expect_true(go_cache_usable(mk_env(), "Homo sapiens", "CC", "ENTREZID"))

    # get_organism() returns NA for an OrgDb it cannot identify (an
    # AnnotationHub OrgDb, say); that used to abort with
    # "missing value where TRUE/FALSE needed"
    expect_false(go_cache_usable(mk_env(), NA_character_, "BP", "ENTREZID"))
    expect_false(go_cache_usable(mk_env(), "Homo sapiens", "BP", NA_character_))

    # nothing cached yet
    expect_false(go_cache_usable(new.env(), "Homo sapiens", "BP", "ENTREZID"))
    expect_false(go_cache_usable(mk_env(goAnno = FALSE), "Homo sapiens", "BP", "ENTREZID"))

    # a different organism or keytype must not be served from the cache
    expect_false(go_cache_usable(mk_env(), "Mus musculus", "BP", "ENTREZID"))
    expect_false(go_cache_usable(mk_env(), "Homo sapiens", "BP", "SYMBOL"))
})

test_that("the GO cache survives an ontology switch", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("GO.db")

    odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    gene <- head(AnnotationDbi::keys(odb, "ENTREZID"), 150)
    run <- function(ont) {
        enrichGO(gene, OrgDb = "org.Hs.eg.db", ont = ont,
                 pvalueCutoff = 1, qvalueCutoff = 1)
    }

    GO_Env <- get_GO_Env()
    ## start from a known cache state -- drop only the annotation cache this
    ## test exercises, so the hybrid GSON cache other tests build survives
    rm(list = intersect(c("goAnno", "keytype", "organism"),
                        ls(GO_Env, all.names = TRUE)),
       envir = GO_Env)

    bp <- run("BP")
    expect_s4_class(bp, "enrichResult")

    ## switching ontology used to evict the annotation and pay a full
    ## mapIds rebuild; it must now be a cache hit with identical results
    expect_identical(as.data.frame(bp), as.data.frame(run("BP")))
    expect_s4_class(mf <- run("MF"), "enrichResult")
    expect_identical(as.data.frame(bp), as.data.frame(run("BP")))
    expect_identical(as.data.frame(mf), as.data.frame(run("MF")))
    expect_s4_class(all <- run("ALL"), "enrichResult")
    expect_identical(as.data.frame(bp), as.data.frame(run("BP")))

    ## the cached annotation is the all-ontology superset
    expect_true(exists("goAnno", envir = GO_Env, inherits = FALSE))
    goAnno <- get("goAnno", envir = GO_Env)
    expect_true(is.factor(goAnno$ONTOLOGYALL))
    expect_true(setequal(unique(as.character(goAnno$ONTOLOGYALL)), c("BP", "CC", "MF")))
})

test_that("enrichGO() can restrict the annotation by evidence code (#160)", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("GO.db")
    skip_if_not_installed("AnnotationDbi")

    odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    gene <- head(AnnotationDbi::keys(odb, "ENTREZID"), 150)
    run <- function(...) {
        enrichGO(gene, OrgDb = "org.Hs.eg.db", ont = "BP",
                 pvalueCutoff = 0.9, qvalueCutoff = 0.9, ...)
    }

    full <- run()
    expect_s4_class(full, "enrichResult")

    # keeping only non-electronic annotations changes the annotation, and so the
    # result; it is not merely a subset of the unfiltered terms
    filt <- run(evidence = c("IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "EXP",
                             "HDA", "HEP", "HMP", "HTP", "IBA", "IC"))
    expect_s4_class(filt, "enrichResult")
    expect_true(nrow(as.data.frame(filt)) > 0)
    expect_false(identical(as.data.frame(full), as.data.frame(filt)))

    # an evidence code the OrgDb does not use leaves nothing, and says so
    expect_error(run(evidence = "NOT_AN_EVIDENCE_CODE"), "no GO annotation left")

    # a filtered run must not poison the cached, unfiltered path
    expect_identical(as.data.frame(full), as.data.frame(run()))
})

test_that("hybrid annotation keeps unmapped input IDs in the result object", {
    skip_if_not_installed("org.Hs.eg.db")

    res <- enrichGO(
        c("TP53", "FAKE_GENE_X"),
        OrgDb = "org.Hs.eg.db",
        keyType = "SYMBOL",
        ont = "BP",
        universe = c("TP53", "BRCA1", "EGFR"),
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 1
    )

    expect_s4_class(res, "enrichResult")
    expect_identical(res@keytype, "SYMBOL")
    expect_true("FAKE_GENE_X" %in% res@gene)
    expect_false(any(grepl("^(E:|N:)", as.character(res@result$geneID))))
})

test_that("hybrid annotation reports an implicit annotated universe", {
    skip_if_not_installed("org.Hs.eg.db")

    res <- expect_message(
        suppressWarnings(enrichGO(
            c("TP53", "FAKE_GENE_X"),
            OrgDb = "org.Hs.eg.db",
            keyType = "SYMBOL",
            ont = "BP",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 1
        )),
        "No universe supplied; using"
    )

    expect_s4_class(res, "enrichResult")
    expect_identical(res@keytype, "SYMBOL")
    expect_true(length(res@universe) > 0)
    expect_false(any(grepl("^(E:|N:)", as.character(res@result$geneID))))
})

## -- hybrid GO data cache ---------------------------------------------------
## build_hybrid_go_data() maps every OrgDb key to Entrez and merges the native
## annotations. With no universe the input is all keys, so the result depends
## only on (organism, ont, keytype) -- but compareCluster() calls enrichGO once
## per cluster and used to rebuild it every time (7.7 s build + 7.2 s restore
## per cluster on human GO BP, PBMC 10k, 2026-10-18). These tests pin that the
## cache is reused for the universe-free path and never poisoned by a
## universe-scoped or evidence-filtered run.

test_that("hybrid_cache_usable requires an exact ontology match", {
    env <- new.env(parent = emptyenv())
    expect_false(clusterProfiler:::hybrid_cache_usable(env, "Homo sapiens", "BP", "SYMBOL"))
    env$hybrid_go <- list()
    env$hybrid_organism <- "Homo sapiens"
    env$hybrid_keytype <- "SYMBOL"
    env$hybrid_ont <- "BP"
    expect_true(clusterProfiler:::hybrid_cache_usable(env, "Homo sapiens", "BP", "SYMBOL"))
    expect_false(clusterProfiler:::hybrid_cache_usable(env, "Homo sapiens", "MF", "SYMBOL"))
    expect_false(clusterProfiler:::hybrid_cache_usable(env, "Homo sapiens", "BP", "ENTREZID"))
    expect_false(clusterProfiler:::hybrid_cache_usable(env, "Mus musculus", "BP", "SYMBOL"))
})

test_that("hybrid path never triggers when every input ID maps", {
    skip_if_not_installed("org.Hs.eg.db")

    env <- clusterProfiler:::get_GO_Env()
    for (what in c("hybrid_go", "hybrid_organism", "hybrid_keytype", "hybrid_ont")) {
        if (exists(what, envir = env, inherits = FALSE)) {
            rm(list = what, envir = env)
        }
    }
    ## All three symbols map, so enrichGO must take the #805 fast path
    ## (map to ENTREZID once, run the Entrez analysis) and must NOT build the
    ## hybrid annotation at all -- no native IDs, no hybrid data, no hybrid
    ## cache write.
    res <- suppressWarnings(suppressMessages(enrichGO(
        c("TP53", "BRCA1", "EGFR"), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
        ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1
    )))
    expect_s4_class(res, "enrichResult")
    expect_false(exists("hybrid_go", envir = env, inherits = FALSE))
    expect_identical(res@keytype, "SYMBOL")
})

test_that("hybrid GO data is built once and reused only for mixed IDs", {
    skip_if_not_installed("org.Hs.eg.db")

    env <- clusterProfiler:::get_GO_Env()
    for (what in c("hybrid_go", "hybrid_organism", "hybrid_keytype", "hybrid_ont")) {
        if (exists(what, envir = env, inherits = FALSE)) {
            rm(list = what, envir = env)
        }
    }
    ## One ID that cannot be mapped to Entrez: this is the only case the hybrid
    ## annotation exists for, and the case where the cache is allowed to act.
    mixed <- c("TP53", "BRCA1", "EGFR", "FAKE_GENE_X")
    run <- function(...) suppressWarnings(suppressMessages(enrichGO(
        mixed, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
        ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, ...
    )))

    first <- run()
    expect_true(exists("hybrid_go", envir = env, inherits = FALSE))
    cached <- get("hybrid_go", envir = env)

    second <- run()
    expect_identical(get("hybrid_go", envir = env), cached)
    expect_equal(as.data.frame(first)$ID, as.data.frame(second)$ID)
    expect_equal(as.data.frame(first)$pvalue, as.data.frame(second)$pvalue)

    ## A universe-scoped run builds from gene + universe: it must not replace
    ## the shared cache with a gene-specific object.
    universe_run <- run(universe = c("TP53", "BRCA1"))
    expect_s4_class(universe_run, "enrichResult")
    expect_identical(get("hybrid_go", envir = env), cached)

    ## An evidence-filtered annotation must not poison the unfiltered cache.
    evidence_run <- run(evidence = "IEA")
    expect_s4_class(evidence_run, "enrichResult")
    expect_identical(get("hybrid_go", envir = env), cached)
})

## -- vectorized gene-set label mapping --------------------------------------
## restore_hybrid_result() and restore_entrez_result_ids() used to map each
## gene set separately; on a BP result that is ~5,900 sets / 390k ID slots and
## measured 7.5 s per call (13 clusters: ~55 s of unname/index work). The
## vectorized helper must be byte-identical, including names and empty sets.
test_that("map_gene_set_labels matches the per-set loop exactly", {
    ids_list <- list(
        g1 = c("a", "b", "c"),
        g2 = character(),
        g3 = c("b", "d"),
        g4 = "e",
        g5 = c("a", "a")
    )
    labels <- stats::setNames(c("A", "B", "C", "D", "E"), c("a", "b", "c", "d", "e"))
    old_sets <- lapply(ids_list, function(ids) unique(unname(labels[as.character(ids)])))
    old_gene_id <- vapply(ids_list, function(ids) {
        paste0(unique(unname(labels[as.character(ids)])), collapse = "/")
    }, character(1))
    got <- clusterProfiler:::map_gene_set_labels(ids_list, labels)
    expect_identical(got$sets, old_sets)
    expect_identical(got$gene_id, old_gene_id)

    ## No gene sets at all must not produce a broken structure.
    empty <- clusterProfiler:::map_gene_set_labels(list(), labels)
    expect_identical(empty$sets, list())
    expect_identical(empty$gene_id, character())
})

test_that("hybrid_ids_for_input keeps every Entrez match and marks unmapped ids", {
    mapping <- data.frame(
        source = c("A", "A", "B", "C"),
        entrez = c("1", "2", "3", "4"),
        stringsAsFactors = FALSE
    )
    got <- clusterProfiler:::hybrid_ids_for_input(c("A", "B", "Z"), mapping, "SYMBOL")
    expect_equal(got, c("E:1", "E:2", "E:3", "N:SYMBOL:Z"))
    ## Duplicated input ids are collapsed, and an empty input stays empty.
    expect_equal(clusterProfiler:::hybrid_ids_for_input(c("A", "A"), mapping, "SYMBOL"),
                 c("E:1", "E:2"))
    expect_identical(clusterProfiler:::hybrid_ids_for_input(character(), mapping, "SYMBOL"),
                     character())
})
