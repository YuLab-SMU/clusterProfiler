## The narrative error-injection suite.
##
## What these tests protect is not the injector -- it is the *claim about the
## safety check*. Each case states an error kind and whether the check is expected
## to see it, so a change that silently widens or narrows that claim fails here.

injection_fixture <- function() {
    b <- structure(
        list(
            regulators = data.frame(
                regulator = c("HIF1A", "MYC", "ESR1", "PGR"),
                source = "tf", term = "t", named = TRUE, level = "B",
                padj_within_source = c(1e-5, 1e-3, 0.5, 0.6),
                setSize = 10L, Count = 5L,
                significant = c(TRUE, TRUE, FALSE, FALSE),
                direction = c("up", "up", "down", "down"),
                stringsAsFactors = FALSE),
            diagnostics = list(counts = data.frame(
                source = "tf", coverage_tested = 0.9, coverage_floor = NA_real_)),
            ## ESR1 and PGR are both unadmitted *and* input genes. That is the
            ## case the check used to miss, so the fixture has to contain it.
            input = list(genes = c("HIF1A", "MYC", "ESR1", "PGR", "GAPDH"),
                         design = list(background = character(0)))),
        class = c("evidenceBundle", "list"))
    structure(
        list(verdicts = data.frame(
                 regulator = c("HIF1A", "MYC", "ESR1", "PGR"),
                 verdict = c("supported", "single_source", "absent", "absent"),
                 max_level = "B", padj = c(1e-5, 1e-3, 0.5, 0.6),
                 direction = c("up", "up", "down", "down"),
                 stringsAsFactors = FALSE),
             evidence = b),
        class = c("adjudication", "list"))
}

coverage_fixture <- function(counts, n_input = 100L, sources = "perturbagen,tfunion") {
    structure(
        list(verdicts = data.frame(
                 regulator = c("HIF1A", "MYC", "ESR1", "PGR"),
                 verdict = c("supported", "single_source", "absent", "absent"),
                 max_level = "B", padj = c(1e-5, 1e-3, 0.5, 0.6),
                 direction = c("up", "up", "down", "down"),
                 stringsAsFactors = FALSE),
             evidence = b),
        class = c("adjudication", "list"))
}

test_that("the compliant narrative names only admitted entities", {
    adj <- injection_fixture()
    text <- .mechanism_narrative(adj)
    expect_type(text, "character")
    expect_gt(nchar(text), 0)
    ## The C6 control is only a control if it really is clean. If this ever
    ## reports something, every false-positive number in the report is wrong.
    expect_length(.mechanism_detect(text, adj), 0)
})

test_that("an unadmitted regulator is caught even when it is an input gene", {
    ## C1. ESR1 is graded absent and is also a differentially expressed gene --
    ## the exact shape that used to pass untouched, because being an input gene
    ## exempted the symbol from the check.
    adj <- injection_fixture()
    text <- .mechanism_corrupt(.mechanism_narrative(adj), "C1", "ESR1")
    expect_true("ESR1" %in% .mechanism_detect(text, adj))
    ## And the same for the other one, so this is not one lucky symbol.
    text2 <- .mechanism_corrupt(.mechanism_narrative(adj), "C1", "PGR")
    expect_true("PGR" %in% .mechanism_detect(text2, adj))
})

test_that("an admitted regulator is not a candidate for injection", {
    ## Injecting an admitted entity would score the check on a case that is
    ## supposed to pass, so the candidate list has to exclude them.
    cand <- .mechanism_injection_candidates(injection_fixture())
    expect_false("HIF1A" %in% cand)
    expect_false("MYC" %in% cand)
    expect_true(all(c("ESR1", "PGR") %in% cand))
})

test_that("a fabricated symbol is caught and is not drawn from the data", {
    ## C3. Generated rather than taken from the bundle: nothing in the evidence
    ## can explain these, so a miss is unambiguously the check's.
    adj <- injection_fixture()
    text <- .mechanism_corrupt(.mechanism_narrative(adj), "C3", "ZQ01")
    expect_true("ZQ01" %in% .mechanism_detect(text, adj))
    expect_false(grepl("ZQ01", paste(adj$evidence$input$genes, collapse = " ")))
})

test_that("symbols of the P53 shape are not mistaken for list markers", {
    ## The caps filter used to drop any single letter followed by digits, on the
    ## reading that such a token is always a panel label. P53, P63, P68, P300 and
    ## F10 are real symbols and were dropped with it, so a narrative could name
    ## p53 as the mechanism and the check would not see it. Found by injecting
    ## every unadmitted regulator of the MCF-7 example: five of nineteen misses.
    adj <- injection_fixture()
    for (sym in c("P53", "P300", "F10")) {
        text <- .mechanism_corrupt(.mechanism_narrative(adj), "C3", sym)
        expect_true(sym %in% .mechanism_detect(text, adj))
    }
    ## The case the rule exists for still holds: a letter and one digit is a
    ## panel label, and reporting it would be noise.
    panel <- .mechanism_corrupt(.mechanism_narrative(adj), "C3", "A4")
    expect_length(.mechanism_detect(panel, adj), 0)
})

test_that("corruption only appends, leaving the original text readable", {
    adj <- injection_fixture()
    base <- .mechanism_narrative(adj)
    for (kind in c("C1", "C3")) {
        out <- .mechanism_corrupt(base, kind, "ZQ07")
        expect_identical(substr(out, 1, nchar(base)), base)
        expect_gt(nchar(out), nchar(base))
    }
    ## C6 is the control and must come back untouched.
    expect_identical(.mechanism_corrupt(base, "C6"), base)
    expect_error(.mechanism_corrupt(base, "C9", "X"), "unknown injection kind")
})

test_that("the report keeps every kind separate and shows the control at zero", {
    adj <- injection_fixture()
    rep <- .mechanism_injection_report(adj)
    expect_s3_class(rep, "data.frame")
    expect_identical(sort(rep$kind), c("C1", "C3", "C6"))
    ## One row per kind, never merged: an average would hide a kind the check
    ## cannot see at all.
    expect_identical(anyDuplicated(rep$kind), 0L)
    ## Both in-scope kinds are caught in full on this fixture.
    expect_equal(rep$rate[rep$kind == "C1"], 1)
    expect_equal(rep$rate[rep$kind == "C3"], 1)
    ## The control reports nothing. Anything else means the checker flags clean
    ## text, which would make every other number here meaningless.
    expect_equal(rep$detected[rep$kind == "C6"], 0L)
    ## C6 injected nothing, so it has no rate -- not a rate of zero.
    expect_true(is.na(rep$rate[rep$kind == "C6"]))
})

test_that("a bundle with nothing admitted yields no narrative and no report", {
    ## Nothing to say is not the same as a clean result, and the report must not
    ## invent an empty narrative to inject into.
    adj <- injection_fixture()
    adj$verdicts$verdict <- "absent"
    expect_null(.mechanism_narrative(adj))
    expect_equal(nrow(.mechanism_injection_report(adj)), 0)
})

## ---- E3: coverage visibility -------------------------------------------------

coverage_bundle <- function(counts, n_input = 100L) {
    b <- structure(
        list(
            regulators = data.frame(regulator = "GR", source = "perturbagen",
                                    term = "D1", named = TRUE, level = "A",
                                    padj_within_source = 1e-23, setSize = 329L,
                                    Count = counts[1], significant = TRUE,
                                    stringsAsFactors = FALSE),
            perturbagen = data.frame(
                regulator = "GR", term = paste0("D", seq_along(counts)),
                method = "shRNA", Count = counts, setSize = 329L,
                padj_within_source = rep(1e-23, length(counts)),
                stringsAsFactors = FALSE),
            provenance = list(upstream = "perturbagen", evidence_sha = "abc"),
            input = list(genes = paste0("g", seq_len(n_input)),
                         design = list(background = character(0)))),
        class = c("evidenceBundle", "list"))
    ## The real bundle carries these; the count line is printed only when it can.
    ## Set here so the test can assert that the count and the coverage line both
    ## appear -- they answer different questions and neither replaces the other.
    attr(b$perturbagen, "queried") <- 1L
    attr(b$perturbagen, "considered") <- 1L
    b
}

coverage_fixture <- function(counts, n_input = 100L, sources = "perturbagen,tfunion") {
    structure(
        list(verdicts = data.frame(
                 regulator = "GR", sources = sources, n_sources = 2L,
                 max_level = "A", verdict = "supported", padj = 1e-23,
                 family = FALSE, evidence_kind = "measured",
                 direction = "up", measured_direction = "mixed",
                 perturbation_effect = NA_character_, writable = "mechanism",
                 stringsAsFactors = FALSE),
             evidence = coverage_bundle(counts, n_input), min_independent = 2,
             coverage_floor = NULL),
        class = "adjudication")
}

test_that("coverage is reported next to the dataset count, not in place of it", {
    ## The count line says how many datasets were queried. It does not say how
    ## much of the caller's list any of them speaks to, which is what a level-A
    ## grade rests on. Airway: GR at A on 23 of 896 genes.
    thin <- capture.output(print(coverage_bundle(counts = 23, n_input = 896)))
    expect_true(any(grepl("coverage:", thin)))
    expect_match(paste(thin, collapse = " "), "median 2\\.6%")

    ## The count line is still there: the two answer different questions.
    expect_true(any(grepl("measured:", thin)))
})

test_that("coverage is shown for every bundle, thin or not", {
    ## This is a description, not a warning. Suppressing it when healthy would
    ## make its absence ambiguous again -- "nothing to report" or "not measured".
    healthy <- capture.output(print(coverage_bundle(counts = c(400, 350), n_input = 1000)))
    expect_match(paste(healthy, collapse = " "), "median 37\\.5%")
    expect_match(paste(healthy, collapse = " "), "none below")
})

test_that("the coverage threshold can be switched off", {
    old <- options(clusterProfiler_coverage_note = 0)
    withr::defer(options(old))
    expect_null(.mechanism_measured_coverage(coverage_bundle(counts = 2)))
})

test_that("a supported entity leaning on thin measured evidence is named", {
    ## Airway GR: supported at level A on 2.6% of the input. The verdict is not
    ## wrong -- the dataset does show the effect -- but it is not the same claim
    ## as one resting on most of the list, and the print should say which it is.
    adj <- coverage_fixture(counts = 23, n_input = 896)
    out <- capture.output(print(adj))
    expect_true(any(grepl("thin:", out)))
    expect_match(paste(out, collapse = " "), "GR \\(2\\.6% of input\\)")
    ## The verdict itself is untouched.
    expect_identical(adj$verdicts$verdict, "supported")
})

test_that("an entity supported only by collections has no coverage to report", {
    ## Reporting "0%" here would be wrong rather than conservative: there is no
    ## measured evidence to cover the list with.
    adj <- coverage_fixture(counts = 23, n_input = 896, sources = "tfunion")
    out <- capture.output(print(adj))
    expect_false(any(grepl("thin:", out)))
})

test_that("a healthy bundle prints no thin line", {
    adj <- coverage_fixture(counts = 700, n_input = 896)
    out <- capture.output(print(adj))
    expect_false(any(grepl("thin:", out)))
})
