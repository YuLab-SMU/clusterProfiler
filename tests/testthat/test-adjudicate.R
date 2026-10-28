## The adjudication object, its print method, and the dispatch in interpret().
##
## Offline: the fixtures are hand-written GMT files, so nothing here downloads.

test_that("adjudicate() returns an adjudication whose columns are the frozen set", {
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("f", 1:5)),
                     log2FoldChange = c(rep(3, 5), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    adj <- adjudicate(ev)

    expect_s3_class(adj, "adjudication")
    expect_setequal(names(adj),
                    c("verdicts", "evidence", "diagnostics", "parameters", "provenance"))

    ## Frozen order. Appending is allowed; renaming or reordering is not,
    ## because every consumer reads this frame by name and the documentation
    ## promises the shape.
    expect_identical(names(adj$verdicts),
                     c("regulator", "sources", "n_sources", "max_level", "verdict",
                       "padj", "family", "evidence_kind", "direction",
                       "measured_direction", "perturbation_effect", "writable"))

    ## The evidence travels with the verdict, so a saved adjudication can be
    ## read next to what produced it.
    expect_s3_class(adj$evidence, "evidenceBundle")
    expect_identical(adj$evidence$regulators, ev$regulators)
})

test_that("the fingerprint identifies the evidence, not the moment", {
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("f", 1:5)),
                     log2FoldChange = c(rep(3, 5), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    a <- adjudicate(ev)
    Sys.sleep(0.1)
    b <- adjudicate(ev)

    ## Two runs over the same evidence must hash alike, or the field is a
    ## timestamp wearing a fingerprint's name.
    expect_identical(a$provenance$evidence_sha, b$provenance$evidence_sha)
    ## ... and the time must still be recorded, just not inside the hash.
    expect_false(identical(a$provenance$time, NA))

    ## A magnitude change that does not cross zero is not different evidence:
    ## no verdict can see it, so the fingerprint must not change either.
    de$log2FoldChange[1] <- 4
    ev_same <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                         design = list(treatment = "none")))
    expect_identical(a$provenance$evidence_sha,
                     adjudicate(ev_same)$provenance$evidence_sha)

    ## Flipping a sign moves a gene between the halves, which the directional
    ## layers do read -- so that is different evidence.
    de$log2FoldChange[1] <- -4
    ev2 <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                     design = list(treatment = "none")))
    expect_false(identical(a$provenance$evidence_sha,
                           adjudicate(ev2)$provenance$evidence_sha))
})

test_that("print() summarises the verdicts and names no per-entity rows", {
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("f", 1:5)),
                     log2FoldChange = c(rep(3, 5), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    adj <- adjudicate(ev)
    txt <- paste(utils::capture.output(print(adj)), collapse = "\n")

    expect_match(txt, "<adjudication>")
    ## All four verdicts are named, including the ones with nothing in them:
    ## "supported 0" is a finding, silence is not.
    for (v in c("supported", "single_source", "absent", "untested")) {
        expect_match(txt, v)
    }
    expect_match(txt, "sha")
})

test_that("an adjudication is not a data.frame, which is why the class exists", {
    ## `interpret()` takes a caller's own data.frame as a documented input type,
    ## so the two have to be distinguishable. A bare frame would collide.
    use_fixture_cache()
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("f", 1:5)),
                     log2FoldChange = c(rep(3, 5), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    ev <- suppressWarnings(evidence(de, universe = de$gene, upstream = "tf",
                                    design = list(treatment = "none")))
    adj <- adjudicate(ev)
    expect_false(is.data.frame(adj))
    expect_true(is.data.frame(adj$verdicts))
    ## consensus_findings() takes either, so callers need not unpack by hand.
    expect_identical(consensus_findings(adj), consensus_findings(adj$verdicts))
})

test_that("adjudicate() refuses anything that is not an evidenceBundle", {
    expect_error(adjudicate(data.frame(regulator = "X")),
                 "must be an evidenceBundle")
    expect_error(adjudicate(list()), "must be an evidenceBundle")
})

test_that("a carried verdict table wins over a fresh adjudication", {
    ## The bug this guards: interpret() re-adjudicated with the default
    ## `min_independent`, so a caller who adjudicated at a stricter setting had
    ## their decision silently overridden -- the one thing that making
    ## adjudication a separate step is supposed to prevent.
    ##
    ## Built by hand rather than through evidence(): the shared fixture's sets
    ## hold five to eight genes, and against those nothing can ever be
    ## significant, so a real bundle here would make the assertion vacuous.
    bundle <- structure(
        list(regulators = data.frame(
                 regulator = c("FOXO3", "FOXO3"), source = c("tf", "perturbagen"),
                 term = c("FOXO3_TARGET_GENES", "D1"),
                 named = c(TRUE, TRUE), level = c("B", "A"),
                 padj_within_source = c(1e-6, 1e-6), setSize = c(8L, 8L),
                 Count = c(5L, 5L), family = c(FALSE, FALSE),
                 evidence_kind = c("measured", "measured"),
                 direction = c("up", "up"), dataset = c(NA, "d"), system = c(NA, "x"),
                 method = c(NA, "knockout"), significant = c(NA, NA),
                 effect_size = c(NA_real_, NA_real_), in_system = c(NA, FALSE),
                 up_rise = c(NA_integer_, 0L), up_fall = c(NA_integer_, 6L),
                 up_sign_padj = c(NA_real_, 1e-6),
                 down_rise = c(NA_integer_, 0L), down_fall = c(NA_integer_, 0L),
                 down_sign_padj = c(NA_real_, NA_real_), stringsAsFactors = FALSE),
             diagnostics = list(counts = data.frame(
                 source = "tf", n_input = 5L, n_tested = 5L, n_sig_covered = 5L,
                 uncovered_sig = 0L, coverage_tested = 1, coverage_floor = NA_real_,
                 status = "graded", stringsAsFactors = FALSE)),
             input = list(genes = paste0("g", 1:5), universe = paste0("g", 1:8),
                          design = list(coverage_floor = NULL, species = "Homo sapiens"),
                          up = paste0("g", 1:5), down = character(0)),
             provenance = list(upstream = "tf", release = "2026.2.Hs")),
        class = c("evidenceBundle", "list"))

    ## At the default the regulator is supported and admitted.
    expect_identical(adjudicate(bundle)$verdicts$verdict, "supported")
    expect_equal(nrow(evidence_admitted(bundle)), 2)

    ## A stricter adjudication -- what a caller gets from min_independent = 3 --
    ## must be what evidence_admitted() uses when it is handed over.
    strict <- adjudicate(bundle, min_independent = 3)$verdicts
    expect_identical(strict$verdict, "single_source")
    expect_equal(nrow(evidence_admitted(bundle, verdicts = strict)), 2)
    ## And a table that admits nothing must admit nothing.
    none <- strict
    none$verdict <- "absent"
    expect_equal(nrow(evidence_admitted(bundle, verdicts = none)), 0)
    ## Both rows of the bundle are excluded, one per source.
    expect_equal(nrow(evidence_excluded(bundle, verdicts = none)), 2)
})
