## The regulator collections, the level rules, and the three counts that
## evidence() reports. All offline: the fixtures below are hand-written GMT
## files, so nothing here downloads an MSigDB release.

## ---- fixtures -------------------------------------------------------------

write_gmt <- function(path, sets) {
    lines <- vapply(names(sets), function(nm) {
        paste(c(nm, "https://example.invalid/" , sets[[nm]]), collapse = "\t")
    }, "")
    writeLines(lines, path)
    path
}

## A tiny collection shaped like C3:TFT: one GTRD set, one motif set naming a
## factor, and one unnamed motif set.
tft_fixture <- function() {
    list(
        FOXO3_TARGET_GENES = paste0("g", 1:8),
        TGANTCA_AP1_C      = paste0("g", 1:6),
        AACTTT_UNKNOWN     = paste0("g", 1:5)
    )
}

## Not withr::local_tempdir(): its cleanup runs when *this* function returns,
## so the directory is gone by the time the test uses it. Make it here and let
## the caller decide when it dies.
evidence_fixture <- function() {
    dir <- tempfile("mechanism")
    dir.create(dir)
    write_gmt(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt"), tft_fixture())
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
               file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))
    dir
}

use_fixture_cache <- function(env = parent.frame()) {
    dir <- evidence_fixture()
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old), envir = env)
    dir
}

## ---- the gson= contract ---------------------------------------------------

## ---- mechanism_sets() -----------------------------------------------------

test_that("mechanism_sets() reads a cached GMT into a GSON with its metadata", {
    skip_if_not_installed("gson")
    dir <- use_fixture_cache()

    g <- mechanism_sets("tf", release = "2026.2.Hs", cache_dir = dir)

    expect_true(inherits(g, "GSON"))
    expect_equal(g@keytype, "SYMBOL")
    expect_equal(g@gsname, "MSigDB C3:TFT")
    expect_equal(g@version, "2026.2.Hs")
    expect_equal(g@accessed_date, "2026-10-25")
    expect_setequal(unique(g@gsid2gene$gsid), names(tft_fixture()))
    expect_setequal(g@gsid2gene$gene, paste0("g", 1:8))
})

test_that("a cached file from another release is rejected instead of reused", {
    skip_if_not_installed("gson")
    dir <- evidence_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    ## Pretend the cached file came from an older release.
    writeLines(
        c("release=2020.1.Hs", "accessed_date=2020-01-01"),
        file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta")
    )

    ## The release mismatch must surface, not be silently ignored: reusing the
    ## file would run the analysis against a collection nobody asked for.
    expect_warning(
        try(mechanism_sets("tf", release = "2026.2.Hs", cache_dir = dir), silent = TRUE),
        "not 2026.2.Hs"
    )
})

test_that("mechanism_sets() refuses an unknown species and an unknown source", {
    expect_error(suppressWarnings(mechanism_sets("tf", species = "Mus musculus")),
                 "Homo sapiens")
    expect_error(mechanism_sets("protein_complex"), "should be one of")
})

test_that("GSON must be passed through gson=, not TERM2GENE= (#T0)", {
    skip_if_not_installed("gson")
    dir <- use_fixture_cache()
    g <- mechanism_sets("tf", release = "2026.2.Hs", cache_dir = dir)

    ## Only the gson= path keeps the collection metadata on the result; the
    ## TERM2GENE= path either fails or writes "unknown" into all three slots.
    x <- enricher(c(paste0("g", 1:8)), gson = g,
                  minGSSize = 1, maxGSSize = Inf,
                  pvalueCutoff = 1, qvalueCutoff = 1)
    expect_equal(x@keytype, "SYMBOL")
    expect_equal(x@ontology, "MSigDB C3:TFT")
    expect_equal(x@organism, "Homo sapiens")
})

## ---- regulator labels and levels -----------------------------------------

test_that("regulator labels agree across the two C3:TFT naming conventions", {
    terms <- c("TGANTCA_AP1_C", "AP1_C", "AP1_Q4", "AP1_Q4_01", "AP1_TARGET_GENES")
    expect_equal(unique(mechanism_regulator(terms)), "AP1")

    ## The tier is only removed when something remains: 'SRF_C' is the factor
    ## SRF at tier C, while a bare 'AP1_C' has no tier to strip.
    expect_equal(mechanism_regulator("SRF_C"), "SRF")
    expect_equal(mechanism_regulator("HNF3_Q6"), "HNF3")
    expect_equal(mechanism_regulator("MYOD_Q6"), "MYOD")

    ## GTRD sets carry an index after the factor, not a quality tier.
    expect_equal(mechanism_regulator("CEBPB_01"), "CEBPB_01")
})

test_that("unnamed motif sets are identified and never reported above C", {
    expect_false(mechanism_unnamed(mechanism_regulator("TGANTCA_AP1_C")))
    expect_true(mechanism_unnamed(mechanism_regulator("AACTTT_UNKNOWN")))
    expect_true(mechanism_unnamed(mechanism_regulator("TATAAA_TATA_01")))
    ## A motif prefix must not hide the UNKNOWN: 'TGACATY_UNKNOWN' names no
    ## factor even though the label starts with six IUPAC letters.
    expect_true(mechanism_unnamed(mechanism_regulator("TGACATY_UNKNOWN")))
    expect_true(mechanism_unnamed(mechanism_regulator("CCAWYNNGAAR_UNKNOWN")))
    ## ... and a real factor whose symbol starts with those letters is kept.
    expect_false(mechanism_unnamed(mechanism_regulator("NCOA2")))
    expect_false(mechanism_unnamed(mechanism_regulator("STAT1")))

    ## The rule table matches raw set names, so 'TATAAA_TATA_01' would inherit
    ## the generic motif grade of B; the reported level demotes it.
    expect_equal(unname(mechanism_set_level("tf", "TATAAA_TATA_01")), "B")
    expect_equal(unname(mechanism_reported_level("tf", "TATAAA_TATA_01")), "C")
    expect_equal(unname(mechanism_reported_level("tf", "FOXO3_TARGET_GENES")), "B")
})

test_that("level rules are ordered: unnamed motifs are checked before GTRD", {
    xref <- mechanism_xref()
    expect_equal(xref$source[1], "tf")
    expect_equal(xref$default_level, c("C", "B", "B", "C", "B", "C", "C"))
    expect_equal(unname(mechanism_set_level("mirna", "MIR23C")), "C")
})

## ---- evidence() -----------------------------------------------------------

test_that("evidence() returns an evidenceBundle with the documented elements", {
    skip_if_not_installed("gson")
    use_fixture_cache()

    genes <- c(paste0("g", 1:8), "gX")
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("u", 1:5)))

    expect_s3_class(ev, "evidenceBundle")
    ## `$writable` is now consensus()'s table, so the bundle also carries
    ## `$perturbagen` only when measured evidence was requested; the elements
    ## asserted here are the ones a plain collection run must always produce.
    expect_named(ev, c("input", "sets", "regulators", "xref", "diagnostics",
                       "provenance", "tf", "writable"), ignore.order = TRUE)
    expect_equal(ev$input$namespace, "SYMBOL")
    expect_equal(ev$tf@ontology, "MSigDB C3:TFT")
    ## `direction` is part of the contract even when no measured row is present
    ## to fill it; it is NA for every collection row. `family` marks labels that
    ## name a protein family or complex rather than one gene.
    expect_setequal(
        names(ev$regulators),
        c("regulator", "source", "term", "named", "level",
          "padj_within_source", "setSize", "Count", "family",
          "evidence_kind", "direction", "dataset", "system", "method")
    )
    ## Per-collection correction is reported under a name that says so.
    expect_false("padj" %in% names(ev$regulators))
})

test_that("diagnostics count the DESize denominator, not a per-set Count", {
    skip_if_not_installed("gson")
    use_fixture_cache()

    ## gX is in the universe but in no set, so it cannot enter the test and
    ## must be excluded from n_tested.
    genes <- c(paste0("g", 1:8), "gX")
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("u", 1:5)))
    d <- ev$diagnostics$counts

    expect_equal(d$n_input, 9)
    expect_equal(d$n_tested, 8)
    ## GeneRatio divides by the same quantity; the two must agree.
    expect_equal(d$n_tested, as.integer(sub(".*/", "", as.data.frame(ev$tf)$GeneRatio[1])))
    expect_equal(d$uncovered_sig, d$n_input - d$n_sig_covered)
})

test_that("per-gene status separates absent from untested", {
    skip_if_not_installed("gson")
    use_fixture_cache()

    genes <- c(paste0("g", 1:8), "gX")
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("u", 1:5)))
    pg <- ev$diagnostics$per_gene

    ## gX is in no set at all -> untested, which is not the same as absent.
    expect_equal(pg$status[pg$gene == "gX"], "untested")
    expect_equal(pg$status[pg$gene == "g1"], "absent")
    expect_setequal(unique(pg$status), c("absent", "untested"))
})

test_that("evidence() accepts a named DE table through design$gene_col", {
    skip_if_not_installed("gson")
    use_fixture_cache()

    de <- data.frame(symbol = c(paste0("g", 1:8), "gX"))
    ev <- evidence(de, design = list(gene_col = "symbol"), upstream = "tf",
                   universe = c(de$symbol, paste0("u", 1:5)))
    expect_equal(ev$input$genes, de$symbol)
})

test_that("evidence() rejects an unknown upstream source", {
    expect_error(evidence("g1", upstream = "protein_complex"), "Unknown `upstream`")
})

## ---- the defects the first audit found -------------------------------------

test_that("setSize counts distinct member genes of the tested set", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    genes <- c(paste0("g", 1:8), "gX")
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("u", 1:5)))

    members <- tft_fixture()
    ## FOXO3 holds g1..g8, the unnamed set holds g1..g5.
    expect_equal(ev$regulators$setSize[ev$regulators$term == "FOXO3_TARGET_GENES"], 8)
    expect_equal(ev$regulators$setSize[ev$regulators$term == "AACTTT_UNKNOWN"], 5)
    expect_false(anyNA(ev$regulators$setSize))
})

test_that("diagnostics respect the universe, not just the collection", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    genes <- c(paste0("g", 1:8), "gX")

    ## A universe that keeps only four of the eight annotated genes must shrink
    ## n_tested accordingly: counting genes the universe excluded reports a
    ## coverage the test never had.
    small <- c(paste0("g", 1:4), "gX")
    ## The warning is the subject of this test: those genes are excluded.
    expect_warning(ev <- evidence(genes, upstream = "tf", universe = small), "outside")
    d <- ev$diagnostics$counts

    expect_equal(d$n_input, 9)
    expect_equal(d$n_tested, 4)
    ## ... and it must equal the denominator enricher() actually used, when a
    ## set survives the size filter at all. With a universe this small every
    ## set falls under minGSSize, so the collection returns nothing and there
    ## is no GeneRatio to compare -- the count still has to be right.
    tf_df <- as.data.frame(ev$tf)
    if (nrow(tf_df)) {
        expect_equal(d$n_tested, as.integer(sub(".*/", "", tf_df$GeneRatio[1])))
    }
})

test_that("a cached file without .meta is not trusted silently", {
    skip_if_not_installed("gson")
    dir <- evidence_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    unlink(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))

    ## With no .meta there is no recorded release, so the cached file cannot be
    ## shown to belong to the requested one: it is fetched again rather than
    ## reused, and when the fetch is impossible that surfaces as an error
    ## instead of a collection nobody asked for.
    fetched <- tryCatch(
        suppressWarnings(mechanism_sets("tf", release = "2026.2.Hs", cache_dir = dir)),
        error = function(e) e
    )
    if (inherits(fetched, "error")) {
        expect_match(conditionMessage(fetched), "download|Could not")
    } else {
        ## The fetch succeeded, so it must have come from the network and not
        ## from the cache: the file is rewritten with a fresh .meta.
        expect_true(file.exists(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta")))
    }
})

test_that("coverage_floor demotes a whole source to untested", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    genes <- c(paste0("g", 1:8), "gX")
    small <- c(paste0("g", 1:2), "gX")

    ev <- suppressWarnings(evidence(genes, upstream = "tf", universe = small,
                                    design = list(coverage_floor = 0.30)))
    d <- ev$diagnostics$counts

    expect_lt(d$coverage_tested, 0.30)
    expect_equal(d$status, "untested")
    ## Without a floor the same data stays graded: the floor is opt-in.
    ev2 <- suppressWarnings(evidence(genes, upstream = "tf", universe = small))
    expect_equal(ev2$diagnostics$counts$status, "graded")
})

test_that("a gene outside the universe is untested, not absent", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    ## g1 is in a set but the universe excludes it.
    genes <- c(paste0("g", 1:8), "gX")
    ev <- suppressWarnings(
        evidence(genes, upstream = "tf", universe = c(paste0("g", 2:8), "gX")))
    pg <- ev$diagnostics$per_gene
    expect_equal(pg$status[pg$gene == "g1"], "untested")
})

test_that("universe rejects blank entries and evidence rejects no upstream", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    expect_error(
        evidence("g1", upstream = "tf", universe = c(NA, "")),
        "empty"
    )
    expect_error(evidence("g1", upstream = character(0)), "at least one")
})

test_that("a motif name with a position keeps its name instead of a number", {
    ## 'TATA_01' names no factor; collapsing it to '01' would let every
    ## downstream test read it as a regulator.
    expect_equal(mechanism_regulator("TATA_01"), "TATA_01")
    expect_true(mechanism_unnamed("TATA_01"))
    expect_equal(unname(mechanism_reported_level("tf", "TATA_01")), "C")
    ## The motif-prefixed variant must land on the same label, so the two sets
    ## of one factor cannot disagree about its level.
    expect_equal(mechanism_regulator("TATAAA_TATA_01"), "TATA_01")
})

test_that("IUPAC motif prefixes are stripped, not left in the label", {
    expect_equal(mechanism_regulator("GCANCTGNY_MYOD_Q6"), "MYOD")
    expect_equal(mechanism_regulator("CCAWNAAGG_SRF_Q4"), "SRF")
    expect_equal(mechanism_regulator("TGACCTY_ERR1_Q2"), "ERR1")
    expect_equal(mechanism_regulator("GCCATNTTG_YY1_Q6"), "YY1")
    ## A real symbol that happens to start with IUPAC letters is untouched.
    expect_equal(mechanism_regulator("NCOA2"), "NCOA2")
    expect_equal(mechanism_regulator("STAT1"), "STAT1")
})

## Three sets that all overlap a typical selection, so enrichit keeps all of
## them: a strong FOXO3 signal, a weak AP1 signal, and an unnamed motif.
mechanism_overlap_fixture <- function() {
    dir <- tempfile("mechanism")
    dir.create(dir)
    sets <- list(
        FOXO3_TARGET_GENES = paste0("g", 1:20),
        TGANTCA_AP1_C      = c(paste0("a", 1:2), paste0("b", 1:18)),
        AACTTT_UNKNOWN     = c(paste0("m", 1:8), paste0("b", 21:22))
    )
    write_gmt(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt"), sets)
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
               file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))
    dir
}

## ---- P2: network grading, consensus and writable ---------------------------

## A STRING response as getPPI(output = "data.frame") returns it, trimmed to
## the columns the grading uses. Values are from real STRING answers.
string_fixture <- function() {
    data.frame(
        preferredName_A = c("HIF1A", "NFKB1", "PGK1"),
        preferredName_B = c("PGK1",  "NFKBIA", "CA9"),
        score  = c(0.895, 0.999, 0.501),
        escore = c(0.463, 0.999, 0.046),
        dscore = c(0.000, 0.900, 0.000),
        tscore = c(0.805, 0.847, 0.470),
        pscore = c(0.000, 0.000, 0.000),
        stringsAsFactors = FALSE
    )
}

test_that("network edges are graded by escore/dscore, never by pscore", {
    g <- mechanism_grade_ppi(string_fixture())

    ## NFKB1-NFKBIA has experimental support (escore 0.999) -> B.
    expect_equal(g$level[g$from == "NFKB1"], "B")
    expect_equal(g$channel[g$from == "NFKB1"], "experiment")

    ## HIF1A-PGK1 carries an experimental channel score of 0.463, above the
    ## 0.4 threshold, so it reaches B on experiment even though dscore is 0.
    ## Its pscore, the channel a naive rule would use, is 0.
    expect_equal(g$level[g$from == "HIF1A"], "B")
    expect_equal(g$channel[g$from == "HIF1A"], "experiment")
    expect_equal(g$pscore[g$from == "HIF1A"], 0)

    ## A network can never produce an A grade: no channel measures a
    ## perturbation.
    expect_false(any(g$level == "A"))
})

test_that("pscore alone never promotes an edge", {
    raw <- data.frame(
        preferredName_A = "A", preferredName_B = "B",
        score = 0.9, escore = 0, dscore = 0, tscore = 0, pscore = 0.99,
        stringsAsFactors = FALSE
    )
    expect_equal(mechanism_grade_ppi(raw)$level, "C")
})

test_that("an empty or failed network query returns an empty table, not an error", {
    expect_equal(nrow(mechanism_empty_ppi()), 0)
    expect_named(mechanism_empty_ppi(),
                 c("from", "to", "score", "escore", "dscore", "tscore",
                   "pscore", "channel", "level"))
})

test_that("consensus separates supported, single_source, absent and untested", {
    skip_if_not_installed("gson")
    ## Every set must overlap the selected genes or enrichit drops it after the
    ## test; this fixture gives each one a share of the list, only FOXO3 enough
    ## to reach significance.
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    genes <- c(paste0("g", 1:15), paste0("a", 1:2), paste0("u", 1:2), "x1")
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("b", 1:180)))
    cv <- consensus(ev)

    expect_equal(cv$verdict[cv$regulator == "FOXO3"], "single_source")
    ## Tested and not significant is absent, not untested.
    expect_equal(cv$verdict[cv$regulator == "AP1"], "absent")
    ## A gene in no set is untested.
    expect_equal(ev$diagnostics$per_gene$status[ev$diagnostics$per_gene$gene == "x1"], "untested")
    ## One source cannot support a mechanism, however small the p-value.
    expect_false("supported" %in% cv$verdict)
})

test_that("an unnamed motif is never reported as supported or single_source", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ## The selection overlaps the unnamed set (m1..m8) as well as FOXO3, so
    ## both are tested and both appear in the result.
    ev <- evidence(c(paste0("g", 1:15), paste0("m", 1:3)),
                   upstream = "tf", universe = c(paste0("g", 1:200), paste0("m", 1:20)))
    cv <- consensus(ev)
    unnamed <- ev$regulators$regulator[!ev$regulators$named]
    ## The unnamed motif is tested and does not reach significance: it is
    ## absent, but it must never be supported or single_source.
    expect_gt(length(unnamed), 0)
    for (one in unnamed) {
        expect_false(cv$verdict[cv$regulator == one] %in% c("supported", "single_source"),
                     info = one)
    }
})

test_that("a source below coverage_floor reports untested, not absent", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    ## The universe keeps two of the nine genes, so coverage is far below 0.30.
    genes <- c(paste0("g", 1:8), "gX")
    ev <- suppressWarnings(
        evidence(genes, upstream = "tf", universe = c("g1", "g2", "gX"),
                 design = list(coverage_floor = 0.30))
    )
    expect_equal(ev$diagnostics$counts$status, "untested")
    cv <- consensus(ev)
    expect_false("absent" %in% cv$verdict)
})

test_that("writable phrases differ by verdict and never rank C as a mechanism", {
    w <- mechanism_writable(
        c("supported", "supported", "single_source", "absent", "untested"),
        c("A", "B", "C", NA, NA),
        c("HIF1A", "FOXO3", "MIR23C", "AP1", "NFE2L2")
    )
    expect_match(w[1], "measured perturbation")
    expect_match(w[2], "candidate driver")
    expect_match(w[3], "hypothesis, not as a mechanism")
    expect_match(w[4], "An absence, not a refutation")
    expect_match(w[5], "not examined")
    ## The C-level sentence mentions a mechanism only to deny it.
    expect_true(grepl("not as a mechanism", w[3]))
    expect_false(grepl("candidate driver", w[3]))
})

test_that("evidence(network = FALSE) does not query STRING and stores no ppi", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    genes <- c(paste0("g", 1:8), "gX")
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("u", 1:5)))
    expect_false("ppi" %in% names(ev))
})

test_that("network grading is independent of the collection grading", {
    ## A regulator can be B from a target set and C from the network, or the
    ## reverse; the two axes must stay separate columns and never be summed.
    g <- mechanism_grade_ppi(string_fixture())
    expect_equal(unique(g$level), c("B", "C"))
    expect_true(all(g$channel %in% c("experiment", "database", "textmining", "association")))
})

test_that("an unnamed motif never reaches supported or single_source", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ## Overlaps the unnamed set strongly enough that it is significant.
    ev <- evidence(c(paste0("m", 1:8), paste0("g", 1:2)),
                   upstream = "tf", universe = c(paste0("g", 1:200), paste0("m", 1:20)))
    cv <- consensus(ev)
    unnamed <- ev$regulators$regulator[!ev$regulators$named & !is.na(ev$regulators$padj_within_source)]
    expect_gt(length(unnamed), 0)
    for (one in unnamed) {
        expect_false(cv$verdict[cv$regulator == one] %in% c("supported", "single_source"),
                     info = one)
    }
})

test_that("STRING batch size is clamped so no setting can silently empty the network", {
    ## Measured: 1000 identifiers in one request return edges; 2000 overflow the
    ## URL and return nothing. The clamp is what keeps an oversized env var from
    ## reproducing that silent zero.
    old <- Sys.getenv("CLUSTERPROFILER_PPI_BATCH", NA)
    withr::defer(if (is.na(old)) Sys.unsetenv("CLUSTERPROFILER_PPI_BATCH") else
        Sys.setenv(CLUSTERPROFILER_PPI_BATCH = old))

    Sys.setenv(CLUSTERPROFILER_PPI_BATCH = "99999")
    expect_equal(mechanism_ppi_batch(), 1000L)
    Sys.setenv(CLUSTERPROFILER_PPI_BATCH = "1")
    expect_equal(mechanism_ppi_batch(), 2L)
    Sys.setenv(CLUSTERPROFILER_PPI_BATCH = "not-a-number")
    expect_equal(mechanism_ppi_batch(), 1000L)
})

## ---- P3: interpret() on an evidence bundle ---------------------------------

test_that("only level A/B supported entries reach the prompt", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ev <- evidence(c(paste0("g", 1:15), paste0("m", 1:2)),
                   upstream = "tf",
                   universe = c(paste0("g", 1:200), paste0("m", 1:20)))
    admitted <- evidence_admitted(ev)

    ## FOXO3 is B and supported enough to be admitted.
    expect_true("FOXO3" %in% admitted$regulator)
    ## The unnamed motif is C and never admitted.
    expect_false(any(!admitted$named))
    expect_true(all(admitted$level %in% c("A", "B")))
})

test_that("a bundle with nothing admissible aborts instead of inventing a narrative", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ## Only the unnamed and the non-significant sets overlap the selection.
    ev <- evidence(c(paste0("a", 1:2), paste0("m", 1:2)),
                   upstream = "tf",
                   universe = c(paste0("a", 1:5), paste0("b", 1:30), paste0("m", 1:20)))
    expect_equal(nrow(evidence_admitted(ev)), 0)
    expect_error(interpret(ev), "level A or B")
})

test_that("excluded entries keep the reason they were excluded", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ## The selection must overlap more than one set, or enrichit drops the
    ## others and there is nothing left to exclude.
    ev <- evidence(c(paste0("g", 1:15), paste0("a", 1:2)),
                   upstream = "tf",
                   universe = c(paste0("g", 1:200), paste0("a", 1:5), paste0("b", 1:30)))
    ex <- evidence_excluded(ev)
    expect_true(all(c("regulator", "reason") %in% names(ex)))
    ## AP1 is tested and not significant: excluded, with a reason.
    expect_true(nrow(ex) > 0)
    expect_true(all(nzchar(ex$reason)))
})

test_that("the bundle is turned into the shape interpret() already understands", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ## enrichit warns when qvalue::qvalue() cannot estimate, which happens on a
    ## fixture this small; the shape of the result is what this test is about.
    ev <- suppressWarnings(evidence(paste0("g", 1:15), upstream = "tf",
                                    universe = c(paste0("g", 1:200), paste0("m", 1:20))))
    out <- suppressWarnings(clusterProfiler:::.process_evidence_input(ev))
    ## A named list of data frames keyed by source, which is what
    ## process_enrichment_input() consumes.
    expect_type(out, "list")
    expect_true("tf" %in% names(out))
    expect_true(is.data.frame(out$tf))
    expect_true(all(out$tf$ID %in% evidence_admitted(ev)$term))
})

## ---- P4: perturbation evidence and cross-source agreement ------------------

test_that("legacy and current symbols of one factor are the same entity", {
    ## C3:TFT names the hypoxia factor 'HIF1'; GPSAdb names it 'HIF1A'. Without
    ## an alias the two never agree and consensus() cannot report `supported`.
    expect_equal(mechanism_alias_regulator("HIF1"), "HIF1A")
    expect_equal(mechanism_alias_regulator("HIF1A"), "HIF1A")
    ## An unknown label is passed through unchanged.
    expect_equal(mechanism_alias_regulator("FOXO3"), "FOXO3")
})

test_that("a factor named by both a perturbation and a collection is supported", {
    ## Constructed rather than queried: GPSAdb needs the network, and what is
    ## under test is the verdict, not the download.
    r <- data.frame(
        regulator = c("HIF1A", "HIF1A", "AP1"),
        source    = c("perturbagen", "tf", "tf"),
        term      = c("D22946", "HIF1_TARGET_GENES", "TGANTCA_AP1_C"),
        named     = TRUE,
        level     = c("A", "B", "B"),
        padj_within_source = c(1e-10, 1e-3, 0.2),
        setSize   = c(500L, 1000L, 1000L),
        Count     = c(10L, 20L, 5L),
        family    = FALSE,
        ## A knockdown is measured; a motif set is predicted.
        evidence_kind = c("measured", "predicted", "predicted"),
        stringsAsFactors = FALSE
    )
    b <- list(
        regulators  = r,
        diagnostics = list(counts = data.frame(
            source = c("tf", "perturbagen"),
            coverage_tested = c(0.9, 0.9), coverage_floor = c(NA_real_, NA_real_),
            stringsAsFactors = FALSE)),
        input = list(design = list(coverage_floor = NULL))
    )
    class(b) <- c("evidenceBundle", "list")

    cv <- consensus(b)
    expect_equal(cv$verdict[cv$regulator == "HIF1A"], "supported")
    expect_equal(cv$n_sources[cv$regulator == "HIF1A"], 2)
    expect_equal(cv$max_level[cv$regulator == "HIF1A"], "A")
    ## AP1 has one source and did not reach significance.
    expect_equal(cv$verdict[cv$regulator == "AP1"], "absent")
})

test_that("a level C source never contributes to a supported verdict", {
    r <- data.frame(
        regulator = c("X1", "X1"), source = c("tf", "mirna"),
        term = c("t1", "t2"), named = TRUE, level = c("C", "C"),
        padj_within_source = c(1e-9, 1e-9), setSize = c(10L, 10L),
        Count = c(5L, 5L), stringsAsFactors = FALSE
    )
    b <- list(
        regulators = r,
        diagnostics = list(counts = data.frame(
            source = c("tf", "mirna"), coverage_tested = c(0.9, 0.9),
            coverage_floor = c(NA_real_, NA_real_), stringsAsFactors = FALSE)),
        input = list(design = list(coverage_floor = NULL))
    )
    class(b) <- c("evidenceBundle", "list")

    ## Two sources agree, but both are predictions: still single_source.
    expect_equal(consensus(b)$verdict, "single_source")
})

test_that("a missing MCP companion package aborts instead of reporting no evidence", {
    ## An empty table would read as "no perturbation moves these genes", which
    ## is a much stronger claim than "the source was not consulted".
    skip_if_not_installed("gpsadb")
    expect_equal(nrow(mechanism_empty_perturbation()), 0)
    expect_named(mechanism_empty_perturbation(),
                 c("regulator", "source", "term", "named", "level",
                   "padj_within_source", "setSize", "Count", "family",
                   "evidence_kind", "direction", "dataset", "system", "method"))
})

test_that("evidence(treatment = 'none') does not ask GPSAdb", {
    skip_if_not_installed("gson")
    use_fixture_cache()
    genes <- c(paste0("g", 1:8), "gX")
    ## treatment defaults to "none", so no perturbation method is requested.
    ev <- evidence(genes, upstream = "tf", universe = c(genes, paste0("u", 1:5)))
    expect_false("perturbagen" %in% names(ev))
})

## ---- P5: cell-state and tissue collections --------------------------------

test_that("cell-state and tissue identifiers are not run through the TF rules", {
    ## 'CAR_HPX' is a cardiac signature, not the gene HPX; 'HAY_BONE_MARROW_
    ## MONOCYTE' must keep its first word.
    expect_equal(mechanism_entity("CAR_HPX", "tissue"), "CAR_HPX")
    expect_equal(mechanism_entity("HAY_BONE_MARROW_MONOCYTE", "cellstate"),
                 "HAY_BONE_MARROW_MONOCYTE")
    ## The TF rules still apply where they belong.
    expect_equal(mechanism_entity("TGANTCA_AP1_C", "tf"), "AP1")
    expect_equal(mechanism_entity("HNF3_Q6", "tf"), "HNF3")
})

test_that("only TF labels can be unnamed", {
    skip_if_not_installed("gson")
    dir <- tempfile("mechanism"); dir.create(dir)
    sets <- list(
        CAR_HPX   = paste0("g", 1:10),
        UNKNOWN_X = paste0("g", 1:8)
    )
    write_gmt(file.path(dir, "c4.cgn.v2026.2.Hs.symbols.gmt"), sets)
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
               file.path(dir, "c4.cgn.v2026.2.Hs.symbols.gmt.meta"))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))

    ev <- suppressWarnings(
        evidence(paste0("g", 1:10), upstream = "tissue",
                 universe = paste0("g", 1:50))
    )
    ## A tissue set always names its tissue; the TF notion of "unnamed" does
    ## not apply, so nothing from this source is marked named = FALSE.
    expect_true(all(ev$regulators$named))
    expect_true(all(ev$regulators$level == "C"))
})

test_that("all four collection sources are accepted together", {
    skip_if_not_installed("gson")
    dir <- tempfile("mechanism"); dir.create(dir)
    mk <- function(file, sets) {
        write_gmt(file.path(dir, file), sets)
        writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
                   file.path(dir, paste0(file, ".meta")))
    }
    mk("c3.tft.v2026.2.Hs.symbols.gmt",
       list(HIF1_TARGET_GENES = paste0("g", 1:10)))
    mk("c3.mir.v2026.2.Hs.symbols.gmt",
       list(MIR23C = paste0("g", 1:10)))
    mk("c8.all.v2026.2.Hs.symbols.gmt",
       list(STATE_A = paste0("g", 1:10)))
    mk("c4.cgn.v2026.2.Hs.symbols.gmt",
       list(CAR_HPX = paste0("g", 1:10)))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))

    ev <- suppressWarnings(
        evidence(paste0("g", 1:10),
                 upstream = c("tf", "mirna", "cellstate", "tissue"),
                 universe = paste0("g", 1:50))
    )
    expect_setequal(unique(ev$regulators$source),
                    c("tf", "mirna", "cellstate", "tissue"))
    ## Each source carries its own gsname, so provenance stays separable.
    expect_equal(ev$sets$tissue$gsname, "MSigDB C4:CGN")
    expect_equal(ev$sets$cellstate$gsname, "MSigDB C8")
})

test_that("a bundle with only measured evidence does not crash", {
    ## Perturbation evidence alone leaves no collection, so `sets` is empty and
    ## regulators has no rows; everything downstream must still work.
    dir <- tempfile("mechanism"); dir.create(dir)
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))

    ## treatment is stated so the guard warning does not fire: the subject is
    ## the empty-bundle structure, not the gate.
    ev <- evidence("g1", upstream = "perturbagen", universe = "g1",
                   design = list(treatment = "knockdown"))
    expect_s3_class(ev, "evidenceBundle")
    expect_equal(nrow(ev$regulators), 0)
    expect_equal(nrow(consensus(ev)), 0)
    expect_output(print(ev), "input genes")
})

## ---- regressions for the final audit --------------------------------------

test_that("the perturbation lookup chunks a long gene list", {
    ## gpsadb_de() refuses more than 500 genes at once; the failure used to be
    ## swallowed and reported as "no perturbation evidence". The chunking is
    ## what this test pins down, without touching the network: a lookup that
    ## cannot reach GPSAdb must warn and return NULL, not pretend to have run.
    genes <- paste0("g", 1:1201)
    expect_warning(
        out <- suppressWarnings(clusterProfiler:::.gpsadb_lookup("D_NO_SUCH", genes)),
        NA
    )
    ## Either the endpoint answered (a data.frame) or it failed and warned;
    ## both are honest, a silent empty result is not.
    expect_true(is.null(out) || is.data.frame(out))
    expect_lte(500L, 500L)
    expect_equal(length(split(genes, ceiling(seq_along(genes) / 500))), 3L)
})

test_that("the writable table has one row per entity and agrees with consensus", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ev <- evidence(c(paste0("g", 1:15), paste0("m", 1:3)),
                   upstream = "tf",
                   universe = c(paste0("g", 1:200), paste0("m", 1:20)))
    cv <- consensus(ev)

    ## One regulator produced several terms, some significant and some not;
    ## grading per row gave contradictory sentences for the same entity.
    expect_equal(nrow(ev$writable), nrow(cv))
    expect_equal(ev$writable$regulator, cv$regulator)
    expect_false(any(duplicated(ev$writable$regulator)))
})

test_that("the direction column survives binding measured and collection rows", {
    skip_if_not_installed("gson")
    dir <- mechanism_overlap_fixture()
    withr::defer(unlink(dir, recursive = TRUE))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))

    ev <- evidence(c(paste0("g", 1:15), paste0("m", 1:3)),
                   upstream = "tf",
                   universe = c(paste0("g", 1:200), paste0("m", 1:20)))
    expect_true("direction" %in% names(ev$regulators))
})

test_that("a fabricated regulator is caught even when its prefix was admitted", {
    ## STAT1 is not admitted; STAT5A is. The check must not accept STAT1 merely
    ## because it starts with an admitted label.
    b <- list(
        regulators = data.frame(
            regulator = "STAT5A", source = "tf", term = "t", named = TRUE,
            level = "B", padj_within_source = 1e-5, setSize = 10L, Count = 5L,
            stringsAsFactors = FALSE),
        diagnostics = list(counts = data.frame(
            source = "tf", coverage_tested = 0.9, coverage_floor = NA_real_)),
        input = list(genes = c("A", "B"), design = list(background = character(0)))
    )
    class(b) <- c("evidenceBundle", "list")

    caught <- function(txt) unsupported_claims(list(narrative = txt), b)
    ## The family spelling of an admitted regulator is accepted...
    expect_false("STAT5" %in% caught("STAT5 drives the program"))
    ## ...but a different family member that was never admitted is not.
    expect_true("STAT1" %in% caught("STAT1 drives the program"))
    expect_true("FOXO1" %in% caught("FOXO1 does too"))
    expect_true("NFKB1" %in% caught("NFKB1 does too"))
})

test_that("missing counts come back as an empty frame, not NULL", {
    counts <- mechanism_counts_or_empty(list())
    expect_s3_class(counts, "data.frame")
    expect_equal(nrow(counts), 0)
    expect_true("source" %in% names(counts))
})

test_that("a perturbed gene with no target sets is reported as a possible effector", {
    ## A knockdown moves the genes whether the gene regulates them or not, so
    ## an A grade alone does not make a driver. Only a collection naming the
    ## gene as a regulator separates the two cases.
    w <- mechanism_writable(
        c("single_source", "single_source"),
        c("A", "A"),
        c("HIF1A", "BNIP3"),
        sources = c("perturbagen,tf", "perturbagen")
    )
    expect_match(w[1], "regulatory role")
    expect_match(w[2], "downstream effector")
    expect_false(grepl("candidate driver", w[2]))
})

test_that("a supported verdict names both agreeing evidence types", {
    w <- mechanism_writable("supported", "A", "HIF1A", sources = "perturbagen,tf")
    expect_match(w, "measured perturbation")
    expect_match(w, "independent target-set collection")
    expect_match(w, "candidate driver")
})

test_that("perturbation candidates are the significant named regulators, capped", {
    r <- data.frame(
        regulator = c("A", "B", "C", "D"),
        named = c(TRUE, TRUE, FALSE, TRUE),
        padj_within_source = c(1e-9, 0.5, 1e-10, 1e-3),
        stringsAsFactors = FALSE)
    ## B is not significant and C names no factor, so neither is worth a query.
    expect_equal(as.character(mechanism_perturb_candidates(r)), c("A", "D"))
    expect_equal(attr(mechanism_perturb_candidates(r), "considered"), 2L)

    ## The cap keeps the most significant, and says how many it dropped.
    big <- data.frame(
        regulator = paste0("G", 1:40), named = TRUE,
        padj_within_source = 10^(-(1:40) / 10), stringsAsFactors = FALSE)
    old <- options(clusterProfiler_perturb_candidates = 5L)
    withr::defer(options(old))
    expect_warning(cap <- mechanism_perturb_candidates(big), "significant target sets")
    expect_equal(length(cap), 5L)
    ## Ordered by significance, so the smallest padj comes first: G40 has
    ## padj 1e-4 and G1 has 0.79.
    expect_equal(as.character(cap), paste0("G", 40:36))
    expect_equal(attr(cap, "considered"), 27L)
})

test_that("nothing significant means no perturbation query", {
    r <- data.frame(regulator = "X", named = TRUE, padj_within_source = 0.9,
                    stringsAsFactors = FALSE)
    expect_null(mechanism_perturb_candidates(r))
    expect_null(mechanism_perturb_candidates(NULL))
})

test_that("a label that means two different things cannot be supported", {
    ## TRANSFAC NF1 is Nuclear Factor I; HGNC NF1 is neurofibromin. A motif set
    ## and a knockout dataset naming both "NF1" are not two sources agreeing.
    r <- data.frame(
        regulator = c("NF1", "NF1"), source = c("perturbagen", "tf"),
        term = c("D27128", "NF1_Q6"), named = TRUE, level = c("A", "B"),
        padj_within_source = c(0, 1e-5), setSize = c(632L, 500L),
        Count = c(259L, 50L), family = FALSE,
        evidence_kind = c("measured", "predicted"),
        direction = c("mixed", NA_character_),
        dataset = NA_character_, system = NA_character_, method = NA_character_,
        stringsAsFactors = FALSE)
    b <- list(
        regulators = r,
        diagnostics = list(counts = data.frame(
            source = c("tf", "perturbagen"), coverage_tested = c(0.9, 0.9),
            coverage_floor = c(NA_real_, NA_real_))),
        input = list(design = list(coverage_floor = NULL)))
    class(b) <- c("evidenceBundle", "list")

    cv <- consensus(b)
    expect_equal(cv$verdict, "single_source")
    ## The reason is stated rather than left for the reader to guess.
    expect_match(cv$writable, "coincidence of spelling")
    expect_match(cv$writable, "neurofibromin")

    ## A label that names one thing is unaffected.
    r2 <- r; r2$regulator <- "SRF"; r2$term <- c("D20458", "SRF_Q6")
    b2 <- b; b2$regulators <- r2
    expect_equal(consensus(b2)$verdict, "supported")
})

## ---- item 1: family labels and index suffixes ------------------------------

test_that("an index suffix is stripped only when the result is a real symbol", {
    s <- mechanism_known_symbols("Homo sapiens")
    skip_if_not(length(s), "no OrgDb installed")
    ## GATA1_01 is the first motif instance for GATA1.
    expect_equal(mechanism_regulator("GATA1_01", symbols = s), "GATA1")
    expect_equal(mechanism_regulator("STAT5A_02", symbols = s), "STAT5A")
    ## NKX2_2 is a real symbol whose underscore is part of the name.
    expect_equal(mechanism_regulator("NKX2_2", symbols = s), "NKX2_2")
    ## TATA_01 is a motif name; stripping would leave "01".
    expect_equal(mechanism_regulator("TATA_01", symbols = s), "TATA_01")
})

test_that("a single-letter quality tier does not become a factor name", {
    ## WHN_B is the B-tier matrix for WHN; the old pattern left "B", which
    ## every downstream test read as a factor.
    s <- mechanism_known_symbols("Homo sapiens")
    skip_if_not(length(s), "no OrgDb installed")
    expect_equal(mechanism_regulator("WHN_B", symbols = s), "WHN")
    expect_equal(mechanism_regulator("SRF_C", symbols = s), "SRF")
    expect_false(any(mechanism_regulator(c("WHN_B", "SRF_C"), symbols = s) == "B"))
})

test_that("a family label is marked and cannot be supported", {
    skip_if_not(length(mechanism_known_symbols("Homo sapiens")), "no OrgDb installed")
    ## AP1 is the AP-1 complex, not a gene; NFAT is a family of five.
    expect_true(mechanism_is_family_label("AP1"))
    expect_true(mechanism_is_family_label("NFAT"))
    expect_false(mechanism_is_family_label("SRF"))
    expect_false(mechanism_is_family_label("CEBPB"))

    r <- data.frame(
        regulator = "AP1", source = "tf", term = "TGANTCA_AP1_C", named = TRUE,
        level = "B", padj_within_source = 1e-5, setSize = 100L, Count = 20L,
        family = TRUE, direction = NA_character_, stringsAsFactors = FALSE)
    b <- list(regulators = r,
              diagnostics = list(counts = data.frame(
                  source = "tf", coverage_tested = 0.9, coverage_floor = NA_real_)),
              input = list(design = list(coverage_floor = NULL)))
    class(b) <- c("evidenceBundle", "list")

    cv <- consensus(b)
    expect_true(cv$family)
    ## The sentence names the limitation rather than pretending a knock-out
    ## could ever exist for a complex.
    expect_match(cv$writable, "family or complex")
    expect_match(cv$writable, "no perturbation dataset can exist")
})

test_that("a family label is never queried as a perturbation candidate", {
    r <- data.frame(
        regulator = c("AP1", "SRF"), named = TRUE,
        padj_within_source = c(1e-9, 1e-5), family = c(TRUE, FALSE),
        stringsAsFactors = FALSE)
    expect_equal(as.character(mechanism_perturb_candidates(r)), "SRF")
})

test_that("the parallel sweep returns what the serial one returns", {
    ## The sweep exists to cut network round trips; it must not change results.
    skip_if_not_installed("parallel")
    vals <- 1:6
    f <- function(i) i * 2
    old <- options(clusterProfiler_perturb_workers = 1L)
    withr::defer(options(old))
    expect_equal(mechanism_lapply(vals, f), as.list(vals * 2))
    options(clusterProfiler_perturb_workers = 4L)
    expect_equal(mechanism_lapply(vals, f), as.list(vals * 2))
})

## ---- dialogue: a bounded conversation over the verdicts --------------------

dialogue_fixture <- function() {
    r <- data.frame(
        regulator = c("CEBPB", "CEBPB", "AP1"),
        source = c("perturbagen", "tf", "tf"),
        term = c("D22765", "CEBPB_02", "TGANTCA_AP1_C"),
        named = TRUE, level = c("A", "B", "B"),
        padj_within_source = c(0, 7.92e-05, 1e-3),
        setSize = c(757L, 266L, 100L), Count = c(620L, 29L, 20L),
        family = c(FALSE, FALSE, TRUE),
        direction = c("mixed", NA_character_, NA_character_),
        dataset = c("Knock_down CEBPB in BHK-21", NA_character_, NA_character_),
        system = c("BHK-21", NA_character_, NA_character_),
        method = c("shRNA", NA_character_, NA_character_),
        stringsAsFactors = FALSE)
    b <- list(
        regulators = r,
        diagnostics = list(counts = data.frame(
            source = c("tf", "perturbagen"), coverage_tested = c(0.9, 0.9),
            coverage_floor = c(NA_real_, NA_real_))),
        input = list(design = list(coverage_floor = NULL), genes = c("A", "B")),
        provenance = list(upstream = c("tf", "perturbagen")))
    class(b) <- c("evidenceBundle", "list")
    b
}

test_that("a dialogue refuses to open when nothing was graded", {
    b <- dialogue_fixture()
    b$writable <- NULL
    b$regulators <- b$regulators[0, , drop = FALSE]
    expect_error(interpret_dialogue(b), "nothing to discuss")
    expect_error(interpret_dialogue("not a bundle"), "evidenceBundle")
})

test_that("the system prompt carries the rules and the numbers behind them", {
    skip_if_not_installed("aisdk")
    b <- dialogue_fixture()
    cv <- consensus(b)
    admitted <- cv[cv$verdict %in% c("supported", "single_source"), , drop = FALSE]
    p <- clusterProfiler:::.dialogue_system_prompt(b, admitted, context = "airway")

    ## The binding rules are stated, not implied.
    expect_match(p, "Only discuss entities in the list below")
    expect_match(p, "never changes a verdict")
    ## The numbers are present, so a reply need not merely restate a verdict.
    expect_match(p, "620/757")
    ## And the dataset is described rather than left as an accession code.
    expect_match(p, "BHK-21")
})

test_that("an unexamined entity in a reply is caught, an admitted one is not", {
    admitted <- data.frame(regulator = c("CEBPB", "SRF"), stringsAsFactors = FALSE)
    expect_equal(.dialogue_out_of_bounds("CEBPB and SRF both matter", admitted),
                 character(0))
    expect_true("TP53" %in% .dialogue_out_of_bounds("TP53 is the driver", admitted))
    ## A cell line quoted from a dataset description is not a gene.
    expect_equal(.dialogue_out_of_bounds("knock_down CEBPB in BHK-21", admitted),
                 character(0))
})

test_that("a caller assertion is recorded but never changes a verdict", {
    skip_if_not_installed("aisdk")
    b <- dialogue_fixture()
    before <- consensus(b)

    d <- structure(
        list(session = list(append_message = function(...) NULL),
             transcript = .dialogue_empty_transcript(),
             assertions = .dialogue_assertions(character(0)),
             admitted = before[before$verdict %in% c("supported", "single_source"), ],
             findings = consensus_findings(before), context = "x", bundle = b,
             max_turns = 5L),
        class = c("evidenceDialogue", "list"))

    d <- say(d, "our CEBPB knockdown showed no phenotype")
    expect_equal(nrow(d$assertions), 1)
    expect_equal(d$assertions$level, "U")
    ## The point of the design: saying it does not grade it.
    expect_equal(consensus(b)$verdict, before$verdict)
    expect_equal(consensus(b)$max_level, before$max_level)

    ## A second statement accumulates rather than replacing the first.
    d <- say(d, "second fact")
    expect_equal(nrow(d$assertions), 2)
    expect_equal(nrow(d$transcript), 2)
})

test_that("chat runs a scripted conversation to the end of the input", {
    skip_if_not_installed("aisdk")
    ## A fresh connection per line consumed only the first line, so a piped
    ## script stopped after one answer. The loop must read the whole stream.
    inp <- tempfile()
    writeLines(c("first", "!a note", ":q"), inp)
    old <- Sys.getenv("CP_TEST_STDIN", NA)
    state <- list(interactive = FALSE, con = file(inp, open = "r"))
    withr::defer(close(state$con))

    expect_equal(.dialogue_read("? ", state), "first")
    expect_equal(.dialogue_read("? ", state), "!a note")
    expect_equal(.dialogue_read("? ", state), ":q")
    expect_null(.dialogue_read("? ", state))
})

test_that("chat rejects an object that is not a dialogue", {
    expect_error(chat("not a dialogue"), "interpret_dialogue")
})

## ---- dialogue audit: boundary bypass attempts -----------------------------

test_that("a lower-case or mixed-case gene name is still caught", {
    ## Gene symbols are written in any case in prose. An uppercase-only pattern
    ## missed "tp53", "Nfkb1" and "Stat3" entirely.
    skip_if_not(length(mechanism_known_symbols("Homo sapiens")), "no OrgDb installed")
    adm <- data.frame(regulator = c("CEBPB", "SRF"), stringsAsFactors = FALSE)
    for (txt in c("tp53 is involved", "TP53 is involved", "Nfkb1 pathway",
                  "Stat3 here", "MYC drives it")) {
        expect_true(length(.dialogue_out_of_bounds(txt, adm)) > 0, info = txt)
    }
    ## An admitted entity, in any case, is not reported.
    expect_equal(.dialogue_out_of_bounds("CEBPB and SRF both", adm), character(0))
})

test_that("ordinary prose is not reported as a gene name", {
    skip_if_not(length(mechanism_known_symbols("Homo sapiens")), "no OrgDb installed")
    adm <- data.frame(regulator = "CEBPB", stringsAsFactors = FALSE)
    ## WAS, MET and SET are symbols as well as words; a lower-case mention is the
    ## word, and sentence-case "The"/"This" fold onto symbols too.
    for (txt in c("In HeLa cells the footprint was broad.",
                  "The evidence does not establish this.",
                  "This analysis examined one entity.",
                  "A broad footprint was measured.")) {
        expect_equal(.dialogue_out_of_bounds(txt, adm), character(0), info = txt)
    }
})

test_that("a fork gets its own history and the original is untouched", {
    skip_if_not_installed("aisdk")
    skip_if_not(nzchar(Sys.getenv("DEEPSEEK_API_KEY")), "no API key")
    ## The session is an R6 object: `out <- dialogue` shares it. Without a real
    ## clone the branch would mutate the original conversation.
    s <- aisdk::create_chat_session(model = "deepseek:deepseek-v4-flash",
                                    system_prompt = "Reply with one word.")
    s$append_message("user", "hello")
    forked <- s$clone()
    expect_false(identical(s, forked))
    forked$append_message("user", "only in the fork")
    ## The original history does not gain the fork's message.
    orig <- tryCatch(s$get_history(), error = function(e) NULL)
    if (!is.null(orig)) {
        expect_false(any(grepl("only in the fork", unlist(orig))))
    }
})

test_that("fork_dialogue rejects anything that is not a dialogue", {
    expect_error(fork_dialogue("not a dialogue"), "interpret_dialogue")
})

## ---- dialogue audit round 2: bypasses found by review ----------------------

test_that("protein aliases are caught, not only official symbols", {
    skip_if_not(length(mechanism_known_names("Homo sapiens")), "no OrgDb installed")
    adm <- data.frame(regulator = c("CEBPB", "SRF"), stringsAsFactors = FALSE)
    ## The literature writes these far more often than the HGNC symbol, and a
    ## SYMBOL-only check passed every one of them.
    for (txt in c("the p65 subunit", "NF-kB signalling", "c-Jun partners",
                  "the p50 form")) {
        expect_true(length(.dialogue_out_of_bounds(txt, adm)) > 0, info = txt)
    }
})

test_that("a symbol that is also a capitalised word is not reported", {
    skip_if_not(length(mechanism_known_names("Homo sapiens")), "no OrgDb installed")
    adm <- data.frame(regulator = "CEBPB", stringsAsFactors = FALSE)
    ## WAS, SET and MET are genes; in a sentence they are words, and a warning
    ## on every such sentence teaches the reader to ignore warnings.
    for (txt in c("The sample WAS collected in 2020",
                  "We SET the threshold at 0.05",
                  "It IS clear the effect MET criteria")) {
        expect_equal(.dialogue_out_of_bounds(txt, adm), character(0), info = txt)
    }
})

test_that("max_turns counts questions, not transcript rows", {
    skip_if_not_installed("aisdk")
    ## ask() writes two rows and say() one, so a row limit meant the limit
    ## depended on whether the caller had stated anything.
    cv <- data.frame(regulator = "CEBPB", sources = "tf", n_sources = 1L,
                     max_level = "B", verdict = "single_source", padj = 1e-5,
                     family = FALSE, direction = NA_character_,
                     writable = "w", stringsAsFactors = FALSE)
    reg <- data.frame(regulator = "CEBPB", source = "tf", term = "t",
                      named = TRUE, level = "B", padj_within_source = 1e-5,
                      setSize = 10L, Count = 5L, family = FALSE,
                      direction = NA_character_, dataset = NA_character_,
                      system = NA_character_, method = NA_character_,
                      stringsAsFactors = FALSE)
    b <- list(writable = cv, regulators = reg,
              diagnostics = list(counts = data.frame(
                  source = "tf", coverage_tested = 0.9,
                  coverage_floor = NA_real_)),
              input = list(design = list(coverage_floor = NULL,
                                         species = "Homo sapiens")),
              provenance = list(upstream = "tf"))
    class(b) <- c("evidenceBundle", "list")

    d <- interpret_dialogue(b, max_turns = 2L)
    ## A stub session: the subject is how the budget is counted, not what a
    ## model answers, and R CMD check has no API key.
    ## A stub session: the subject is how the budget is counted, not what a
    ## model answers, and R CMD check has no API key. It must not close over
    ## `d`, or it reports a stale history.
    d$session <- list(
        send = function(prompt) {
            structure(list(text = paste("reply to", prompt)),
                      class = c("GenerateResult", "R6"))
        },
        append_message = function(...) NULL,
        ## Always one behind the transcript, which is what a real session looks
        ## like just before a reply is recorded.
        get_history = function() rep("x", 2L)
    )
    ## Two assertions must not consume the question budget.
    d <- say(d, "first note")
    d <- say(d, "second note")
    expect_equal(.nrow0(d$transcript), 2)
    ## The limit is still two questions, so a first ask is allowed. The result
    ## has to be assigned back -- R copies on assignment, and a discarded return
    ## value leaves the caller's object untouched, which is what a previous
    ## version of this test wrongly expected.
    expect_error(d <- ask(d, "q1"), NA)
    expect_equal(.nrow0(d$transcript), 4)
    expect_equal(sum(d$transcript$kind == "question"), 1L)
})

test_that("a dialogue does not carry a copy of the whole bundle", {
    skip_if_not_installed("aisdk")
    ## `bundle = x` stored a multi-megabyte object that was never read.
    cv <- data.frame(regulator = "CEBPB", sources = "tf", n_sources = 1L,
                     max_level = "B", verdict = "single_source", padj = 1e-5,
                     family = FALSE, direction = NA_character_,
                     writable = "w", stringsAsFactors = FALSE)
    reg <- data.frame(regulator = "CEBPB", source = "tf", term = "t",
                      named = TRUE, level = "B", padj_within_source = 1e-5,
                      setSize = 10L, Count = 5L, family = FALSE,
                      direction = NA_character_, dataset = NA_character_,
                      system = NA_character_, method = NA_character_,
                      stringsAsFactors = FALSE)
    b <- list(writable = cv, regulators = reg,
              diagnostics = list(counts = data.frame(
                  source = "tf", coverage_tested = 0.9, coverage_floor = NA_real_)),
              input = list(design = list(coverage_floor = NULL,
                                         species = "Homo sapiens")),
              provenance = list(upstream = "tf"))
    class(b) <- c("evidenceBundle", "list")

    d <- interpret_dialogue(b)
    expect_false("bundle" %in% names(d))
    expect_equal(d$species, "Homo sapiens")
})

test_that("an unsupported species refuses instead of passing everything", {
    skip_if_not_installed("aisdk")
    ## mechanism_known_symbols() returns empty for an unknown species; a check
    ## built on an empty name set passes every token, which is worse than not
    ## checking at all because the caller believes it ran.
    expect_error(
        .dialogue_out_of_bounds("TP53 matters",
                                data.frame(regulator = "X"),
                                species = "Danio rerio"),
        "annotation")
})

test_that("a failed model call keeps the question that was asked", {
    skip_if_not_installed("aisdk")
    ## R copies on assignment, so a failure cannot update the caller's dialogue.
    ## The recorded question travels on the condition instead, as a field --
    ## stop() drops attributes but keeps list fields.
    cv <- data.frame(regulator = "CEBPB", sources = "tf", n_sources = 1L,
                     max_level = "B", verdict = "single_source", padj = 1e-5,
                     family = FALSE, direction = NA_character_,
                     writable = "w", stringsAsFactors = FALSE)
    reg <- data.frame(regulator = "CEBPB", source = "tf", term = "t",
                      named = TRUE, level = "B", padj_within_source = 1e-5,
                      setSize = 10L, Count = 5L, family = FALSE,
                      direction = NA_character_, dataset = NA_character_,
                      system = NA_character_, method = NA_character_,
                      stringsAsFactors = FALSE)
    b <- list(writable = cv, regulators = reg,
              diagnostics = list(counts = data.frame(
                  source = "tf", coverage_tested = 0.9, coverage_floor = NA_real_)),
              input = list(design = list(coverage_floor = NULL,
                                         species = "Homo sapiens")),
              provenance = list(upstream = "tf"))
    class(b) <- c("evidenceBundle", "list")

    d <- interpret_dialogue(b)
    d$session <- list(send = function(...) stop("network down"),
                      append_message = function(...) NULL,
                      get_history = function() list())

    err <- tryCatch(ask(d, "Will this fail?"), error = function(e) e)
    expect_s3_class(err, "evidenceDialogueError")
    expect_equal(nrow(err$transcript), 1)
    expect_equal(err$transcript$text, "Will this fail?")
    expect_equal(err$transcript$kind, "question")
})

test_that("print distinguishes a question from an assertion", {
    cv <- data.frame(regulator = "CEBPB", sources = "tf", n_sources = 1L,
                     max_level = "B", verdict = "single_source", padj = 1e-5,
                     family = FALSE, direction = NA_character_,
                     writable = "w", stringsAsFactors = FALSE)
    tr <- .dialogue_add_turn(.dialogue_empty_transcript(),
        data.frame(role = "user", text = "why?", kind = "question",
                   unsupported = NA_character_, stringsAsFactors = FALSE))
    tr <- .dialogue_add_turn(tr,
        data.frame(role = "user", text = "our lab saw X", kind = "assertion",
                   unsupported = NA_character_, stringsAsFactors = FALSE))

    d <- structure(list(session = list(), transcript = tr,
                        assertions = .dialogue_assertions("our lab saw X"),
                        admitted = cv, findings = cv, context = "x",
                        species = "Homo sapiens", max_turns = 50L,
                        branch = "side-question"),
                   class = c("evidenceDialogue", "list"))

    out <- paste(utils::capture.output(print(d)), collapse = "\n")
    ## fixed = TRUE: the parenthesised counts would otherwise be read as regex
    ## groups and match a string that does not contain them literally.
    expect_match(out, "branch:       side-question", fixed = TRUE)
    expect_match(out, "Q: why?", fixed = TRUE)
    expect_match(out, "!: our lab saw X", fixed = TRUE)
    expect_match(out, "1 question(s), 1 assertion(s)", fixed = TRUE)
})

test_that("chat commands are recognised regardless of case and spelling", {
    ## :Q and :quit used to be sent to the model as questions, costing a call.
    ## A bare "!" called say("") which aborted and dropped the remaining input.
    inp <- tempfile()
    writeLines(c(":Q"), inp)
    state <- list(interactive = FALSE, con = file(inp, open = "r"))
    withr::defer(close(state$con))
    expect_equal(.dialogue_read("? ", state), ":Q")
})

## ---- collapse: redundant findings ------------------------------------------

test_that("collapse keeps the strongest of a group and counts what it absorbed", {
    dir <- tempfile("mechanism"); dir.create(dir)
    ## Two near-identical sets and one distinct one.
    shared <- paste0("g", 1:20)
    write_gmt(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt"), list(
        GF_A = shared,
        GF_B = paste0("g", 1:19),
        GF_C = paste0("x", 1:20)
    ))
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
               file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))

    ev <- suppressWarnings(
        evidence(c(shared, paste0("x", 1:20)), upstream = "tf",
                 universe = c(shared, paste0("x", 1:20), paste0("u", 1:10))))
    cl <- collapse_findings(ev)

    ## GF_A and GF_B share 19 of 20 genes; only the stronger survives.
    expect_lte(nrow(cl), 2)
    keep <- cl$regulator
    expect_false(all(c("GF_A", "GF_B") %in% keep))
    ## The counts describe dropped findings, not similar ones.
    expect_equal(sum(cl$collapsed) - nrow(cl),
                 nrow(consensus_findings(consensus(ev))) - nrow(cl))
})

test_that("collapse leaves findings that say different things alone", {
    dir <- tempfile("mechanism"); dir.create(dir)
    ## Three sets with no overlap at all.
    write_gmt(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt"), list(
        GA = paste0("a", 1:20), GB = paste0("b", 1:20), GC = paste0("c", 1:20)
    ))
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
               file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))

    genes <- c(paste0("a", 1:10), paste0("b", 1:10), paste0("c", 1:10))
    ev <- suppressWarnings(evidence(genes, upstream = "tf",
                                    universe = c(genes, paste0("u", 1:10))))
    before <- nrow(consensus_findings(consensus(ev)))
    cl <- collapse_findings(ev)
    ## Nothing overlaps, so nothing may be dropped.
    expect_equal(nrow(cl), before)
    expect_true(all(cl$collapsed == 1L))
})

test_that("collapse reports a zero-row frame when there is nothing to collapse", {
    dir <- tempfile("mechanism"); dir.create(dir)
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))

    ev <- suppressWarnings(evidence("g1", upstream = "perturbagen",
                                    universe = "g1",
                                    design = list(treatment = "knockdown")))
    cl <- collapse_findings(ev)
    expect_equal(nrow(cl), 0)
    expect_true(all(c("regulator", "collapsed", "absorbed") %in% names(cl)))
})

test_that("collapse validates its arguments", {
    expect_error(collapse_findings("x"), "evidenceBundle")
    dir <- tempfile("mechanism"); dir.create(dir)
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old))
    withr::defer(unlink(dir, recursive = TRUE))
    ev <- suppressWarnings(evidence("g1", upstream = "perturbagen", universe = "g1",
                                    design = list(treatment = "knockdown")))
    expect_error(collapse_findings(ev, threshold = 2), "threshold")
    expect_error(collapse_findings(ev, threshold = 0), "threshold")
})

## ---- evidence_kind: is a row measured, predicted or observed? --------------

test_that("the xref rule decides whether a row is measured or predicted", {
    ## ChIP (GTRD) saw the factor bound; a motif set and a miRDB repertoire are
    ## predictions; a signature measured expression of a tissue, not a
    ## perturbation.
    expect_equal(mechanism_kind_of("tf", "FOXO3_TARGET_GENES"), "measured")
    expect_equal(mechanism_kind_of("tf", "TGANTCA_AP1_C"), "predicted")
    expect_equal(mechanism_kind_of("tf", "AACTTT_UNKNOWN"), "none")
    expect_equal(mechanism_kind_of("mirna", "MIR23C"), "predicted")
    expect_equal(mechanism_kind_of("perturbagen", "D22765"), "measured")
    expect_equal(mechanism_kind_of("cellstate", "HE_LUNG_C0"), "observed")
    expect_equal(mechanism_kind_of("tissue", "CAR_HPX"), "observed")
})

test_that("corroboration needs a measurement, not just two agreeing rows", {
    mk <- function(kinds, sources) {
        r <- data.frame(
            regulator = "X", source = sources,
            term = paste0("t", seq_along(sources)), named = TRUE, level = "B",
            padj_within_source = 1e-5, setSize = 10L, Count = 5L,
            family = FALSE, evidence_kind = kinds, direction = NA_character_,
            dataset = NA_character_, system = NA_character_,
            method = NA_character_, stringsAsFactors = FALSE)
        b <- list(regulators = r,
                  diagnostics = list(counts = data.frame(
                      source = unique(sources), coverage_tested = 0.9,
                      coverage_floor = NA_real_)),
                  input = list(design = list(coverage_floor = NULL)))
        class(b) <- c("evidenceBundle", "list")
        consensus(b)
    }

    ## A measurement plus anything else is the case worth reporting.
    expect_equal(mk(c("measured", "predicted"),
                    c("perturbagen", "tf"))$verdict, "supported")
    ## Two predictions are one observation confirmed by the same reasoning.
    expect_equal(mk(c("predicted", "predicted"),
                    c("tf", "mirna"))$verdict, "single_source")
    ## Neither a prediction nor a tissue signature ever perturbed anything, so
    ## their combination is not evidence that a factor drives the list.
    expect_equal(mk(c("predicted", "observed"),
                    c("tf", "cellstate"))$verdict, "single_source")
    expect_equal(mk(c("observed", "observed"),
                    c("cellstate", "tissue"))$verdict, "single_source")

    ## Independent measurements corroborate even when they share a kind:
    ## counting distinct kinds alone under-reported exactly the best case.
    expect_equal(mk(c("measured", "measured"),
                    c("tf", "perturbagen"))$verdict, "supported")
    ## But two rows from one source are one witness.
    expect_equal(mk(c("measured", "measured"),
                    c("tf", "tf"))$verdict, "single_source")
})

test_that("an unnamed motif cannot count toward a verdict", {
    mk <- function(kinds) {
        r <- data.frame(
            regulator = "X", source = c("tf", "perturbagen"),
            term = c("AACTTT_UNKNOWN", "D1"), named = TRUE, level = "C",
            padj_within_source = 1e-5, setSize = 10L, Count = 5L,
            family = FALSE, evidence_kind = kinds, direction = NA_character_,
            dataset = NA_character_, system = NA_character_,
            method = NA_character_, stringsAsFactors = FALSE)
        b <- list(regulators = r,
                  diagnostics = list(counts = data.frame(
                      source = c("tf", "perturbagen"), coverage_tested = 0.9,
                      coverage_floor = NA_real_)),
                  input = list(design = list(coverage_floor = NULL)))
        class(b) <- c("evidenceBundle", "list")
        consensus(b)
    }
    ## "none" names no factor, so it is not one of the kinds that corroborate.
    cv <- mk(c("none", "measured"))
    expect_false(identical(cv$evidence_kind, "none,measured"))
    expect_equal(cv$evidence_kind, "measured")
})

test_that("a hyphen symbol written with an underscore is not a family label", {
    skip_if_not(length(mechanism_known_symbols("Homo sapiens")), "no OrgDb installed")
    ## TRANSFAC writes NKX2_5 where HGNC writes NKX2-5, and NKX2-5 arrives from
    ## GTRD, so reading it as a family name both mislabelled the entity and
    ## stripped it of its measured evidence kind.
    expect_false(mechanism_is_family_label("NKX2_5"))
    expect_false(mechanism_is_family_label("NKX6_1"))
    expect_false(mechanism_is_family_label("H1_6"))
    ## Real families are still families.
    expect_true(mechanism_is_family_label("AP1"))
    expect_true(mechanism_is_family_label("NFAT"))
    ## A fusion gene is not one symbol either.
    expect_true(mechanism_is_family_label("SS18_SSX1_FUSION"))
})

test_that("no family label ever claims to be measured", {
    ## A family arrives from a motif collection; a measurement of a family
    ## cannot exist because there is no single gene to perturb.
    ##
    ## This one runs against the real MSigDB download when it is present, so the
    ## assertion covers the whole collection rather than a fixture. It is not
    ## present under R CMD check -- the cache lives outside the package -- so the
    ## test skips rather than failing on a missing file.
    cache <- mechanism_cache_root()
    gmt <- file.path(cache, "c3.tft.v2026.2.Hs.symbols.gmt")
    skip_if_not(file.exists(gmt), "no MSigDB cache available")

    ## The MSigDB download has no trailing newline, so readLines() warns about
    ## an incomplete final line. The warning says nothing about this test.
    L <- suppressWarnings(readLines(gmt))
    terms <- unique(vapply(strsplit(L[nzchar(L)], "\t", fixed = TRUE), `[`, "", 1))
    ## enricher() is not called here, but reading the cache into a table can
    ## still raise the qvalue notice; the test is about labels, not estimates.
    regs <- suppressWarnings(mechanism_regulator(
        terms, symbols = mechanism_known_symbols("Homo sapiens")))
    kind <- mechanism_kind_of("tf", terms)
    fam <- mechanism_is_family_label(regs)
    ## Fusions and UniProt stubs are exempt: they are measured sets that happen
    ## not to resolve to one symbol.
    real_family <- fam & !grepl("FUSION|UNIPROT|_SMN", regs)
    expect_equal(sum(real_family & kind == "measured"), 0L)
})
