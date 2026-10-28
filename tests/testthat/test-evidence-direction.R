## A collection big enough for the pathway layer to have any power.
##
## The shared fixture's sets hold five to eight genes, and against it a five-gene
## pathway yields p = 1 however the arithmetic is arranged -- there is no room
## for an over-representation to show. This fixture writes its own collection to
## a temporary cache, so evidence() and the built result read the same genes:
## 500 background genes, a 100-gene pathway, and a 20-gene regulator set that
## shares ten of them.
pathway_fixture <- function(env = parent.frame()) {
    dir <- tempfile("pathway")
    dir.create(dir)
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    write_gmt(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt"),
              list(PATH1 = paste0("u", 1:100),
                   TFX_TARGET_GENES = paste0("u", 1:20)))
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-28"),
               file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old), envir = env)
    list(universe = paste0("u", 1:500),
         pathgenes = paste0("u", 1:100))
}

## An enrichment result over that pathway, built with the same call a caller
## would use.
pathway_result <- function(fx) {
    enricher(fx$pathgenes, universe = fx$universe, gson = mechanism_sets("tf"),
             pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
}

test_that("a DE table carries its direction through to the up/down split", {
    use_fixture_cache()
    ## Five of the eight set genes rise, the rest sit still. The set must be
    ## reported against the up half, not against the list as a whole.
    de <- data.frame(
        gene = c(paste0("g", 1:5), paste0("x", 1:5)),
        log2FoldChange = c(3, 3, 3, 3, 3, -0.2, -0.2, -0.2, -0.2, -0.2),
        stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    expect_true(ev$input$directed)
    r <- ev$regulators
    expect_identical(r$direction[r$term == "FOXO3_TARGET_GENES"][1], "up")
})

test_that("a table with no fold-change column says so instead of running blind", {
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("x", 1:5)),
                     stringsAsFactors = FALSE)
    expect_warning(
        ev <- evidence(de, universe = de$gene, upstream = "tf",
                       design = list(treatment = "none")),
        "no direction")
    expect_false(ev$input$directed)

    ## Naming a column that does not exist is a different thing: the caller
    ## asked for a specific column and it is missing, so guessing again on their
    ## behalf would hide the mistake.
    expect_error(
        evidence(de, universe = de$gene, upstream = "tf",
                 design = list(gene_col = "gene", logfc_col = "logFC")),
        "no column \"logFC\"")
})

test_that("an enrichResult is refused as genes, with the replacement named", {
    use_fixture_cache()
    dir <- getOption("clusterProfiler_cache_dir")
    g <- mechanism_sets("tf")
    ## minGSSize is explicit: enricher()'s default of 10 would drop every
    ## fixture set (they hold five to eight genes) and return NULL, which is an
    ## empty answer rather than an error.
    res <- enricher(paste0("g", 1:5), universe = paste0("g", 1:8), gson = g,
                    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 2)
    expect_s4_class(res, "enrichResult")
    ## The loss is the reason: @result is truncated at pvalueCutoff, so reading
    ## genes back out of it silently drops some. The message has to say what to
    ## do instead, or the caller just looks for another way to pass the result.
    expect_error(evidence(res), "no longer accepts an `enrichResult`")
    expect_error(evidence(res), "pathway = result")
})

test_that("pathway = FALSE leaves the layer alone", {
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("x", 1:5)),
                     log2FoldChange = c(3, 3, 3, 3, 3, -1, -1, -1, -1, -1),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    expect_null(ev$pathways)
})

test_that("identical gene sets collapse to one pathway, aliases kept", {
    ## GO's near-synonyms are the reason: three names, one set of genes. Testing
    ## each name would run the same test three times and report three findings
    ## that are one finding.
    df <- data.frame(
        ID = c("GO:1", "GO:2", "GO:3"),
        Description = c("response to hypoxia", "response to decreased oxygen levels",
                        "response to oxygen levels"),
        geneID = c("A/B/C", "C/A/B", "B/C/A"),
        p.adjust = c(1e-13, 1e-12, 1e-11),
        stringsAsFactors = FALSE)
    out <- clusterProfiler:::.mechanism_pathways_dedupe(df)
    expect_identical(nrow(out), 1L)
    ## The best-corrected name represents the set.
    expect_identical(out$Description[1], "response to hypoxia")
    expect_match(out$aliases[1], "decreased oxygen")
})

test_that("a reused enrichment result supplies pathways without being recomputed", {
    fx <- pathway_fixture()
    pres <- pathway_result(fx)
    expect_s4_class(pres, "enrichResult")
    ## Symbols, because enricher() wrote them; that is the shape the layer wants.
    expect_true(grepl("^u[0-9]+/", as.data.frame(pres)$geneID[1]))

    de <- data.frame(gene = c(paste0("u", 1:10), paste0("d", 1:5)),
                     log2FoldChange = c(rep(3, 10), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = fx$universe, upstream = "tf",
                                    pathway = pres,
                                    design = list(treatment = "none")))
    expect_s3_class(ev$pathways, "data.frame")
    expect_true(nrow(ev$pathways) > 0)
    ## The pathway row says which pathway it is about, not only which regulator.
    expect_true(all(c("pathway", "pathway_description", "regulator",
                      "padj_within_pathway") %in% names(ev$pathways)))
    expect_true("TFX" %in% ev$pathways$regulator)
})

test_that("the pathway layer never becomes a witness in the verdict table", {
    fx <- pathway_fixture()
    pres <- pathway_result(fx)
    ## The TF sets are the same collection the collection layer uses, so the
    ## pathway layer is a re-reading of one source, not a second one.
    de <- data.frame(gene = c(paste0("u", 1:10), paste0("d", 1:5)),
                     log2FoldChange = c(rep(3, 10), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    with_pw <- suppressWarnings(evidence(de, universe = fx$universe, upstream = "tf",
                                         pathway = pres,
                                         design = list(treatment = "none")))
    without <- suppressWarnings(evidence(de, universe = fx$universe, upstream = "tf",
                                         design = list(treatment = "none")))
    ## Assert the layer really produced rows, or the comparisons below are
    ## vacuous.
    expect_true(nrow(with_pw$pathways) > 0)
    ## Derived from the same collections, so counting it as a source would let
    ## one collection corroborate itself -- the failure that made
    ## evidence_biotrigger() raise `supported` without adding a regulator.
    expect_false("pathway" %in% with_pw$regulators$source)
    expect_setequal(unique(with_pw$regulators$source), c("tf"))
    a <- consensus_findings(consensus(with_pw))
    b <- consensus_findings(consensus(without))
    expect_identical(a$n_sources, b$n_sources)
    expect_identical(a$verdict, b$verdict)
})

## ---- L2: the perturbation's direction --------------------------------------

test_that("the sign test reads a side only when the test supports it", {
    ## The old rule required every moved gene to agree. A knockdown moves some
    ## genes up and some down, so on any real list that never held and every
    ## dataset came back `mixed` -- a property of the ruler, not of the data.
    moved <- data.frame(gene = paste0("g", 1:30),
                        log2fc = c(rep(-2, 12), rep(2, 2), rep(0.5, 16)),
                        padj = 1e-6)
    up <- paste0("g", 1:14)
    half <- clusterProfiler:::.mechanism_sign_test(moved, up)
    expect_identical(half$fall, 12L)
    expect_identical(half$rise, 2L)
    expect_lt(half$p, 0.05)
    half$padj <- half$p
    expect_identical(clusterProfiler:::.mechanism_half_direction(half), "down")

    ## Two observations cannot produce p < 0.05, so a half that small must report
    ## no direction rather than a test that always says `mixed`.
    tiny <- clusterProfiler:::.mechanism_sign_test(moved, c("g1", "g15"))
    expect_true(is.na(tiny$p))
    expect_true(is.na(tiny$direction))
})

test_that("perturbation_effect names the mode of action, not agreement", {
    ## Both of these are mechanisms. A boolean "concordant" flag would call the
    ## second one a failure, which it is not: it says the factor restrains these
    ## genes. On the hypoxia contrast that was 7 of 25 supported findings.
    row <- function(direction, method, measured, down = c(0L, 0L, NA_real_)) {
        data.frame(source = c("tf", "perturbagen"),
                   direction = c(direction, measured),
                   method = c(NA, method),
                   down_rise = c(0L, down[1]), down_fall = c(0L, down[2]),
                   down_sign_padj = c(NA, down[3]),
                   stringsAsFactors = FALSE)
    }
    act <- clusterProfiler:::.consensus_perturbation_effect(row("up", "knockout", "down"))
    rep <- clusterProfiler:::.consensus_perturbation_effect(row("up", "knockout", "up"))
    expect_identical(act, "activator")
    expect_identical(rep, "repressor")

    ## Overexpression inverts the expectation: the factor is supplied, so its
    ## targets should rise for an activator.
    expect_identical(
        clusterProfiler:::.consensus_perturbation_effect(row("up", "overexpression", "up")),
        "activator")

    ## No direction, no claim.
    expect_true(is.na(clusterProfiler:::.consensus_perturbation_effect(
        row("up", "knockout", "mixed"))))
    expect_true(is.na(clusterProfiler:::.consensus_perturbation_effect(
        row("up", "unknown_method", "down"))))
})

test_that("a repressor report does not change the verdict", {
    ## The flag is a report, not a condition. Requiring it would re-decide every
    ## verdict the package has already made, and those were not wrong.
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("f", 1:5)),
                     log2FoldChange = c(rep(3, 5), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    cv <- consensus(ev)
    ## Every regulator in this fixture is a collection-only finding, so there is
    ## nothing to report and nothing that could have flipped a verdict.
    expect_true(all(is.na(cv$perturbation_effect)))
    expect_true("perturbation_effect" %in% names(cv))
})
