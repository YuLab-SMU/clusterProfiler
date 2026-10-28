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
