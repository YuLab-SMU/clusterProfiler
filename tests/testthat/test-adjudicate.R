## The adjudication object, its print method, and the dispatch in interpret().
##
## Offline: the fixtures are hand-written GMT files, so nothing here downloads.

adjudication_contract_bundle <- function() {
    r <- data.frame(regulator = "FOXO3", source = c("tf", "perturbagen"),
                    term = c("FOXO3_TARGETS", "D1"), named = TRUE,
                    level = c("B", "A"), padj_within_source = 1e-6,
                    setSize = 8L, Count = 5L, family = FALSE,
                    evidence_kind = c("predicted", "measured"),
                    direction = c("up", "down"), method = c(NA, "knockout"),
                    stringsAsFactors = FALSE)
    structure(list(regulators = r, perturbagen = r[2, , drop = FALSE],
                   diagnostics = list(counts = data.frame(
                       source = c("tf", "perturbagen"), coverage_tested = c(0.4, 1),
                       coverage_floor = NA_real_, status = "graded")),
                   input = list(genes = paste0("g", 1:5), universe = paste0("g", 1:8),
                                up = paste0("g", 1:5), down = character(),
                                design = list(coverage_floor = NULL, species = "Homo sapiens")),
                   provenance = list(upstream = c("tf", "perturbagen"), release = "fixture")),
              class = c("evidenceBundle", "list"))
}

test_that("explicit coverage floors control verdicts and admitted source rows", {
    ev <- adjudication_contract_bundle()
    expect_identical(adjudicate(ev)$verdicts$verdict, "supported")
    legacy <- ev
    legacy$diagnostics$counts$status <- NULL
    expect_identical(adjudicate(legacy)$verdicts$verdict, "supported")
    strict <- adjudicate(ev, coverage_floor = 0.8)
    expect_identical(strict$verdicts$verdict, "single_source")
    expect_identical(strict$verdicts$sources, "perturbagen")
    expect_identical(strict$diagnostics$counts$status, c("untested", "graded"))
    expect_identical(evidence_admitted(ev, verdicts = strict$verdicts)$source, "perturbagen")
    excluded <- evidence_excluded(ev, verdicts = strict$verdicts)
    expect_identical(excluded$source, "tf")
    expect_identical(excluded$reason, "source excluded by adjudication")
    ev$input$design$coverage_floor <- 0.8
    ev$diagnostics$counts$coverage_floor <- 0.8
    expect_identical(adjudicate(ev)$verdicts$verdict, "single_source")
    expect_identical(adjudicate(ev, coverage_floor = 0)$verdicts$verdict, "supported")
    ev$diagnostics$counts$coverage_tested <- 0.1
    expect_identical(adjudicate(ev, coverage_floor = 0.8)$verdicts$verdict, "untested")
})

test_that("fingerprints include actual evidence and decision parameters but not row order", {
    ev <- adjudication_contract_bundle()
    sha <- adjudicate(ev)$provenance$evidence_sha
    changed <- ev
    changed$regulators$padj_within_source <- 0.9
    expect_identical(adjudicate(changed)$verdicts$verdict, "absent")
    expect_false(identical(sha, adjudicate(changed)$provenance$evidence_sha))
    for (field in c("direction", "method", "level", "term")) {
        changed <- ev
        changed$regulators[[field]][1] <- "changed"
        expect_false(identical(sha, .adjudication_sha(changed)), info = field)
    }
    expect_false(identical(sha, adjudicate(ev, min_independent = 3)$provenance$evidence_sha))
    expect_false(identical(sha, adjudicate(ev, coverage_floor = 0.8)$provenance$evidence_sha))
    changed <- ev
    changed$regulators <- changed$regulators[2:1, rev(names(changed$regulators))]
    changed$diagnostics$counts <- changed$diagnostics$counts[2:1, ]
    changed$provenance$called_at <- "a different moment"
    expect_identical(sha, adjudicate(changed)$provenance$evidence_sha)
})

test_that("two views of D1 plus a TF collection are only two witnesses", {
    ev <- adjudication_contract_bundle()
    row <- ev$regulators[2, , drop = FALSE]
    row$source <- "biotrigger"
    ev$regulators <- rbind(ev$regulators, row)
    expect_identical(adjudicate(ev)$verdicts$verdict, "supported")
    expect_identical(adjudicate(ev, min_independent = 3)$verdicts$verdict, "single_source")
})

test_that("adjudication authority retains identity direction and writable strength", {
    cv <- adjudicate(adjudication_contract_bundle(), min_independent = 3)$verdicts
    rows <- jsonlite::fromJSON(.adjudication_prompt_rows(cv))
    fields <- c("regulator", "verdict", "direction", "measured_direction",
                "perturbation_effect", "writable")
    expect_equal(rows[, fields], cv[, fields])
    expect_identical(rows$verdict, "single_source")
    rejected <- cv
    rejected$regulator <- "MYC"
    rejected$verdict <- "absent"
    expect_equal(jsonlite::fromJSON(.adjudication_prompt_rows(rbind(cv, rejected))), rows)
    mixed <- cv
    mixed$direction <- "mixed"
    mixed$measured_direction <- NA_character_
    serialized <- jsonlite::fromJSON(.adjudication_prompt_rows(mixed))
    expect_identical(serialized$direction, "mixed")
    expect_true(is.na(serialized$measured_direction))
    prompt <- .adjudication_declaration(cv)
    expect_match(prompt, "never promote")
    expect_match(prompt, "preserve exact regulator identity")
    expect_match(prompt, "No verdict proves an upstream causal driver")
    expect_match(prompt, "ADJUDICATION_ROWS_JSON")
    expect_match(prompt, cv$writable, fixed = TRUE)
    expect_match(prompt, "not a guarantee of model compliance", fixed = TRUE)
})

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

test_that("pathway attributions reach the prompt as context and only via decided entities", {
    ## The layer re-reads collections the verdicts already used, so it must add
    ## no entity and strengthen no claim. Both guards are checked: an entity the
    ## adjudication did not admit never appears, and the block carries its own
    ## instruction -- putting it in the system prompt would change the bytes of
    ## the four legacy input types.
    pw <- data.frame(
        pathway = c("GO:1", "GO:1", "GO:2"),
        pathway_description = c("response to hypoxia", "response to hypoxia", "glycolysis"),
        pathway_direction = c("up", "up", "up"),
        regulator = c("HIF1A", "MYC", "NOT_ADMITTED"),
        padj_within_pathway = c(1e-5, 1e-3, 1e-4),
        pathway_padj = c(1e-13, 1e-13, 1e-2),
        stringsAsFactors = FALSE)
    adj <- structure(
        list(verdicts = data.frame(regulator = c("HIF1A", "MYC"),
                                   verdict = c("supported", "single_source"),
                                   stringsAsFactors = FALSE),
             evidence = list(pathways = pw)),
        class = c("adjudication", "list"))

    ctx <- .mechanism_pathway_context(adj)
    txt <- paste(ctx, collapse = "\n")
    expect_match(txt, "HIF1A")
    expect_match(txt, "context, not evidence")
    expect_false(grepl("NOT_ADMITTED", txt))
    ## The pathway whose only attribution was unadmitted disappears entirely
    ## rather than appearing with an empty list.
    expect_false(grepl("glycolysis", txt))

    ## The most significant pathways come first, not the alphabetically first.
    ## The layer groups its rows by pathway id, so row order is no guide.
    shuffled <- pw
    shuffled$pathway <- c("GO:0000001", "GO:0000001", "GO:9999999")
    shuffled$pathway_description <- c("worst", "worst", "best")
    shuffled$padj_within_pathway <- c(0.1, 0.1, 1e-20)
    ## The pathway's own enrichment p is what orders the section, so give the
    ## one with the weaker attribution the stronger pathway.
    shuffled$pathway_padj <- c(1e-20, 1e-20, 0.5)
    adj2 <- structure(
        list(verdicts = data.frame(regulator = "HIF1A", verdict = "supported",
                                   stringsAsFactors = FALSE),
             evidence = list(pathways = shuffled)),
        class = c("adjudication", "list"))
    ordered <- .mechanism_pathway_context(adj2)
    expect_lt(regexpr("best", paste(ordered, collapse = " ")),
              regexpr("worst", paste(ordered, collapse = " ")))

    ## Rendered in its own section, next to the terms it is a reading of.
    prompt <- .build_data_sections("TERM1 (p=1e-3)", pathway_context = ctx)
    expect_match(prompt, "Regulator-Pathway Attribution")
    expect_lt(regexpr("Regulator-Pathway Attribution", prompt),
              regexpr("HIF1A", prompt))

    ## Nothing admitted, nothing to say.
    empty <- structure(
        list(verdicts = data.frame(regulator = "X", verdict = "absent",
                                   stringsAsFactors = FALSE),
             evidence = list(pathways = pw)),
        class = c("adjudication", "list"))
    expect_null(.mechanism_pathway_context(empty))
})

test_that("network evidence is reference material and changes no verdict", {
    ## Locked so that nobody later assumes it contributes. Measured on 54
    ## datasets whose true regulator is known: the edge set is real and on target,
    ## and it moved no answer, because no code turns an edge into corroboration.
    use_fixture_cache()
    universe <- c(paste0("g", 1:8), paste0("f", 1:200))
    de <- data.frame(gene = c(paste0("g", 1:5), paste0("f", 1:5)),
                     log2FoldChange = c(rep(3, 5), rep(-1, 5)),
                     stringsAsFactors = FALSE)
    without <- suppressWarnings(evidence(de, universe = universe, upstream = "tf",
                                         design = list(treatment = "none")))
    ## Exercise the network branch without depending on STRING availability.
    local_mocked_bindings(evidence_ppi = function(...) mechanism_empty_ppi(),
                          .package = "clusterProfiler")
    with_net <- suppressWarnings(evidence(de, universe = universe, upstream = "tf",
                                         network = TRUE,
                                         design = list(treatment = "none")))
    expect_identical(without$regulators, with_net$regulators)
    expect_identical(adjudicate(without)$verdicts, adjudicate(with_net)$verdicts)
})

test_that("the adjudication declaration is appended and never substituted", {
    ## The four legacy input types keep the prompt they had byte for byte, so
    ## the declaration may only ever be added on top -- and only for a run that
    ## came through adjudicate().
    base <- .interpretation_system_prompt()
    expect_false(grepl("ADJUDICATION", base))
    expect_identical(.interpretation_system_prompt_for(), base)
    expect_identical(.interpretation_system_prompt_for(character(0)), base)
    expect_identical(.interpretation_system_prompt_for(NULL), base)

    with_decl <- .interpretation_system_prompt_for(.adjudication_declaration())
    ## Appended, not replaced: the base has to survive inside the result.
    expect_identical(substr(with_decl, 1, nchar(base)), base)
    expect_gt(nchar(with_decl), nchar(base))
    ## Every verdict is named, because the point is to stop each one being
    ## upgraded in prose.
    for (v in c("supported", "single_source", "absent", "untested")) {
        expect_match(with_decl, v)
    }
    ## Annotation prompts are built elsewhere and never carry the declaration.
    expect_false(grepl("ADJUDICATION", .annotation_system_prompt("c1")))
})

test_that("measured rows are read independently of their order", {
    ## Two measured datasets for one regulator used to be resolved by taking the
    ## first row: reversing the bundle flipped perturbation_effect between
    ## activator and repressor while the evidence fingerprint stayed identical.
    ev <- adjudication_contract_bundle()
    p <- ev$regulators[2, , drop = FALSE]
    b <- p
    b$source <- "biotrigger"; b$term <- "D2"; b$method <- "overexpression"
    bundle <- ev
    bundle$regulators <- rbind(ev$regulators[1, , drop = FALSE], p, b)
    fwd <- adjudicate(bundle)
    bundle$regulators <- bundle$regulators[3:1, , drop = FALSE]
    rev <- adjudicate(bundle)
    expect_identical(fwd$provenance$evidence_sha, rev$provenance$evidence_sha)
    expect_identical(fwd$verdicts, rev$verdicts)
    ## knockout+down is an activator, overexpression+down is a repressor: the
    ## measured evidence disagrees, so no single effect is claimed.
    expect_true(is.na(fwd$verdicts$perturbation_effect))

    expect_identical(fwd$verdicts$measured_direction, "down")

    ## Opposite measurements agree after normalizing by perturbation method.
    agree <- bundle
    agree$regulators <- rbind(ev$regulators[1, , drop = FALSE], p, b)
    agree$regulators$direction[3] <- "up"
    both <- adjudicate(agree)
    expect_identical(both$verdicts$perturbation_effect, "activator")
    expect_identical(both$verdicts$measured_direction, "mixed")
    agree$regulators <- agree$regulators[3:1, , drop = FALSE]
    expect_identical(adjudicate(agree)$verdicts, both$verdicts)
})

test_that("mixed and unknown measurements do not invent opposing effects", {
    cases <- list(
        list(direction = "mixed", method = "knockout", measured = "mixed"),
        list(direction = NA_character_, method = "knockout", measured = "down"),
        list(direction = "up", method = "unknown_method", measured = "mixed"))
    for (case in cases) {
        ev <- adjudication_contract_bundle()
        other <- ev$regulators[2, , drop = FALSE]
        other$term <- "D2"
        other$direction <- case$direction
        other$method <- case$method
        ev$regulators <- rbind(ev$regulators, other)
        fwd <- adjudicate(ev)
        expect_identical(fwd$verdicts$measured_direction, case$measured)
        expect_identical(fwd$verdicts$perturbation_effect, "activator")
        expect_identical(fwd$verdicts$verdict, "supported")
        ev$regulators <- ev$regulators[3:1, , drop = FALSE]
        rev <- adjudicate(ev)
        expect_identical(rev$verdicts, fwd$verdicts)
        expect_identical(rev$provenance$evidence_sha, fwd$provenance$evidence_sha)
    }
    ev <- adjudication_contract_bundle()
    ev$regulators$direction[2] <- "mixed"
    expect_identical(adjudicate(ev)$verdicts$measured_direction, "mixed")
    expect_true(is.na(adjudicate(ev)$verdicts$perturbation_effect))
})

test_that("mixed collections do not select a directional programme half", {
    for (direction in c("up", "down")) {
        ev <- adjudication_contract_bundle()
        ev$regulators$direction[1] <- direction
        ev$regulators$down_rise <- 0L
        ev$regulators$down_fall <- 12L
        ev$regulators$down_sign_padj <- 1e-6
        other <- ev$regulators[1, , drop = FALSE]
        other$term <- "OTHER_TARGETS"
        other$direction <- "mixed"
        ev$regulators <- rbind(ev$regulators, other)
        fwd <- adjudicate(ev)
        expect_identical(fwd$verdicts$direction, "mixed")
        expect_true(is.na(fwd$verdicts$measured_direction))
        expect_true(is.na(fwd$verdicts$perturbation_effect))
        expect_identical(fwd$verdicts$verdict, "supported")
        ev$regulators <- ev$regulators[3:1, , drop = FALSE]
        expect_identical(adjudicate(ev)$verdicts, fwd$verdicts)
    }
})

test_that("missing or nonfinite down-half statistics leave adjudication unknown", {
    base <- adjudication_contract_bundle()
    base$regulators$direction[1] <- "down"
    base$regulators$down_rise <- 0L
    base$regulators$down_fall <- 12L
    base$regulators$down_sign_padj <- 1e-6
    valid <- adjudicate(base)$verdicts
    expect_identical(valid$measured_direction, "down")
    expect_identical(valid$perturbation_effect, "activator")
    for (field in c("down_rise", "down_fall", "down_sign_padj")) {
        for (value in list(NULL, NA_real_, Inf, c(1, 2))) {
            ev <- base
            ev$regulators[[field]] <- if (length(value) > 1L) {
                rep(list(value), nrow(ev$regulators))
            } else value
            fwd <- adjudicate(ev)
            expect_true(is.na(fwd$verdicts$measured_direction), info = field)
            expect_true(is.na(fwd$verdicts$perturbation_effect), info = field)
            expect_identical(fwd$verdicts$verdict, "supported")
            ev$regulators <- ev$regulators[2:1, , drop = FALSE]
            expect_identical(adjudicate(ev)$verdicts, fwd$verdicts)
        }
    }
})
