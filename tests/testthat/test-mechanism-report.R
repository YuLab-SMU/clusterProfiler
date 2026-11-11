## 判决书层的离线测试。
##
## 这一层最容易出的错不是崩溃，而是**动作与原因不对应**：一条 concept 侧的失败
## 被写成"再补一份数据"，用户就会去做一件不会改变结论的事。所以测试重点放在
## "每个 reason 都有对应动作、且三侧不串"上，而不是只测能不能打印。

.make_pe <- function(level, verdict, reason, driver = NA_character_, pathway = "p") {
    data.frame(pathway = pathway, direction = "up", n_members = 20L, support = 0.8,
               driver = driver, dataset = "D1", cell_line = "A", es = -0.5,
               padj = 1e-3, effect = "activator", n_tests_run = 2L,
               n_cell_lines_supported = 1L, level = level, verdict = verdict,
               reason = reason, statement = "s", stringsAsFactors = FALSE)
}

test_that("every reason pathway_evidence can emit has an action and a side", {
    ## 穷举而不是抽样：`pathway_evidence()` 能产出的 reason 值，每一个都必须能在
    ## 动作表里查到。查不到的会静默变成 NA 动作，读者看到一行"没有下一步"却不知道
    ## 是没有还是漏了。
    reasons <- c("single-cell-line", "direction-conflict", "no-dataset", "too-small",
                 "not-significant", "member-derived", "no-members", "low-support",
                 "no-direction", "no-fold-changes", "no-test-completed", "lookup-failed",
                 "no-perturbation-test", "coherent-refuted")
    tab <- clusterProfiler:::.mechanism_remedies()
    expect_equal(anyDuplicated(tab$reason), 0L)
    expect_true(all(tab$side %in% c("data", "analysis", "concept")))
    for (r in reasons) {
        hit <- clusterProfiler:::.mechanism_remedy_for(r, "unsupported")
        expect_false(is.na(hit$side))
        expect_false(is.na(hit$action))
    }
    ## 等级 A 没有 reason，也要有一句下一步（在独立扰动里验证）。
    hit <- clusterProfiler:::.mechanism_remedy_for(NA_character_, "perturbation-supported")
    expect_false(is.na(hit$side))
})

test_that("no concept-side failure is turned into an experiment", {
    tab <- clusterProfiler:::.mechanism_remedies()
    concept <- tab[tab$side == "concept", , drop = FALSE]
    expect_gt(nrow(concept), 0L)
    expect_true(all(!concept$unlocks))
})

test_that("verdict classes split supported, conflicting and undetermined", {
    pe <- rbind(.make_pe("A", "perturbation-supported", NA),
                .make_pe("B", "direction-conflict", "direction-conflict"),
                .make_pe("C", "unsupported", "too-small"))
    expect_equal(clusterProfiler:::.mechanism_verdict_class(pe),
                 c("supported", "conflicting", "undetermined"))
})

test_that("a single cell line is reported as a gap, with the tissue when known", {
    ## 库里只有 `NR3C1` 这一个写法；用户写 `GR`。别名必须能在两边对上，
    ## 否则报告会对 `GR` 说"库里没有这个驱动者的数据"——一句假的事实陈述。
    inv <- data.frame(perturbed_gene = "NR3C1", dataset_id = "D1",
                      cell_line = "HEK293T", tissue = "Kidney",
                      n_samples = 4, stringsAsFactors = FALSE)
    g <- clusterProfiler:::.mechanism_gap_table(inv, "GR", list(tissue = "Airway"))
    expect_equal(g$n_datasets[1], 1L)
    expect_equal(g$n_cell_lines[1], 1L)
    expect_match(g$gap[1], "no available dataset from Airway")
    expect_match(g$gap[1], "Kidney")
})

test_that("an inventory without perturbed_gene is an error, not 'no data'", {
    bad <- data.frame(driver = "HIF1A", dataset_id = "D1", cell_line = "HeLa",
                      tissue = "Cervix", n_samples = 4, stringsAsFactors = FALSE)
    expect_error(clusterProfiler:::.mechanism_inventory("HIF1A", bad), "perturbed_gene")
})

test_that("the report ranks only fixable actions and keeps concept rows out", {
    pe <- rbind(.make_pe("C", "unsupported", "single-cell-line", pathway = "a"),
                .make_pe("C", "unsupported", "single-cell-line", pathway = "b"),
                .make_pe("C", "coherent-refuted", "not-significant", pathway = "c"))
    rep <- mechanism_report(pe, drivers = character(0))
    expect_s3_class(rep, "mechanismReport")
    expect_equal(rep$next_steps$unlocks[1], 2L)
    expect_true(all(rep$next_steps$side %in% c("data", "analysis")))
    expect_equal(sum(rep$verdicts$class == "undetermined"), 3L)
    expect_equal(sum(rep$verdicts$class == "conflicting"), 0L)
})

test_that("the report does not mistake nonsignificance for exclusion", {
    pe <- .make_pe("C", "coherent-refuted", "not-significant")
    rep <- mechanism_report(pe, drivers = character(0))
    expect_equal(nrow(rep$next_steps), 0L)
    expect_equal(rep$verdicts$class, "undetermined")
    expect_match(rep$summary[length(rep$summary)], "does not rule out")
    expect_false(any(grepl("No experiment would change", rep$summary)))
})

test_that("a blocked preflight stops the report from recommending experiments", {
    pe <- rbind(.make_pe("C", "unsupported", "single-cell-line", pathway = "a"),
                .make_pe("C", "unsupported", "too-small", pathway = "b"))
    pf <- structure(list(tier = "description",
                         reasons = c("few-up-genes: 23 genes have log2FC > 1; below 50, pathways are too small to carry a driver claim."),
                         metrics = data.frame(n_pathways = 2L, n_up = 23L,
                                              median_members = 2, n_cell_lines = 1L),
                         counts = list(n_pathways = 2L, n_up = 23L)),
                    class = "mechanismPreflight")
    rep <- mechanism_report(pe, drivers = character(0), preflight = pf)
    expect_equal(nrow(rep$next_steps), 0L)
    expect_match(paste(rep$summary, collapse = " "), "attribution is not supported")
    expect_match(paste(rep$summary, collapse = " "), "present evidence limit")
    expect_false(any(grepl("No experiment would change", rep$summary)))
    rep2 <- mechanism_report(pe, drivers = character(0))
    expect_gt(nrow(rep2$next_steps), 0L)
})

test_that("print covers evidence sections and the inventory", {
    pe <- rbind(.make_pe("C", "unsupported", "single-cell-line", driver = "GR", pathway = "a"),
                .make_pe("C", "coherent-refuted", "not-significant", pathway = "b"))
    rep <- mechanism_report(pe, inventory = data.frame(
        perturbed_gene = "GR", dataset_id = "D1", cell_line = "HEK293T",
        tissue = "Kidney", n_samples = 3, stringsAsFactors = FALSE))
    out <- capture.output(print(rep))
    expect_true(any(grepl("mechanism report", out)))
    expect_true(any(grepl("^  conflict:", out)))
    expect_true(any(grepl("undetermined", out)))
    expect_true(any(grepl("driver inventory", out)))
    expect_true(any(grepl("next experiment", out)))
})

test_that("a bare too-small verdict still gets an action", {
    ## 键必须落在补救表里：原来写 "too-small-decomp"（全仓只此一处），
    ## match 落空后那一行会静默没有下一步。
    hit <- clusterProfiler:::.mechanism_remedy_for(NA_character_, "too-small")
    expect_false(is.na(hit$side))
    expect_false(is.na(hit$action))
})

test_that("direction-conflict is not sold as a data gap", {
    ## 它是证据反对（读作压制），不能承诺"再补一份数据就能救"。
    hit <- clusterProfiler:::.mechanism_remedy_for("direction-conflict", "direction-conflict")
    expect_equal(hit$side, "concept")
    expect_false(isTRUE(hit$unlocks))
})

test_that("drivers whose target sets overlap are put in one module", {
    ## Unsigned overlap is a descriptive neighbourhood, not indistinguishability.
    net <- data.frame(
        source = c(rep("HIF1A", 5), rep("EPAS1", 4), rep("ATF4", 5)),
        target = c("a","b","c","d","e", "a","b","c","f", "g","h","i","j","k"),
        stringsAsFactors = FALSE)
    m <- clusterProfiler:::.mechanism_module_map(c("HIF1A","EPAS1","ATF4"), net, 0.10)
    expect_setequal(m$HIF1A, c("HIF1A","EPAS1"))
    expect_setequal(m$EPAS1, c("EPAS1","HIF1A"))
    expect_equal(m$ATF4, "ATF4")          # 靶集不重叠，自己一组
    ## 阈值调高就分开
    m2 <- clusterProfiler:::.mechanism_module_map(c("HIF1A","EPAS1","ATF4"), net, 0.9)
    expect_equal(m2$HIF1A, "HIF1A")
})

test_that("the module columns are filled only when a network is supplied", {
    x <- data.frame(pathway = c("p1","p2"), level = c("B","B"),
                    verdict = c("single-dataset","single-dataset"),
                    reason = NA_character_, driver = c("HIF1A","EPAS1"),
                    es = c(-1,-1), padj = c(0.01,0.01), stringsAsFactors = FALSE)
    net <- data.frame(source = c(rep("HIF1A",5), rep("EPAS1",4)),
                      target = c("a","b","c","d","e","a","b","c","f"),
                      stringsAsFactors = FALSE)
    inv <- data.frame(perturbed_gene = c("HIF1A","EPAS1"),
                      cell_line = c("A","B"), tissue = c("T","T"),
                      dataset_id = c("d1","d2"), stringsAsFactors = FALSE)
    none <- mechanism_report(x, inventory = inv)
    expect_true(all(is.na(none$verdicts$module)))
    with_net <- mechanism_report(x, inventory = inv, network = net)
    expect_equal(with_net$verdicts$module[1], "HIF1A, EPAS1")
    expect_equal(with_net$verdicts$module_size[1], 2L)
})

test_that("a broad knockdown is flagged without changing the class", {
    x <- data.frame(pathway = c("p1","p2"), level = c("B","B"),
                    verdict = c("single-dataset","single-dataset"),
                    reason = NA_character_, driver = c("ATF4","XBP1"),
                    es = c(-1,-1), padj = c(0.01,0.01), stringsAsFactors = FALSE)
    inv <- data.frame(perturbed_gene = c("ATF4","XBP1"),
                      n_up = c(2000, 600), n_down = c(2000, 500),
                      cell_line = c("A","B"), tissue = c("T","T"),
                      dataset_id = c("d1","d2"), stringsAsFactors = FALSE)
    r <- mechanism_report(x, inventory = inv, global_effect_de = 3000L)
    expect_match(r$verdicts$global_effect_reason[1], "global-effect")
    expect_true(is.na(r$verdicts$global_effect_reason[2]))
    ## 判定档位不受影响
    expect_equal(r$verdicts$class, c("undetermined","undetermined"))
    expect_true(any(grepl("global_effect_reason", r$summary)))
})

test_that("module and global-effect never change the verdict class", {
    x <- data.frame(pathway = "p1", level = "A",
                    verdict = "perturbation-supported", reason = NA_character_,
                    driver = "HIF1A", es = -1, padj = 0.01, stringsAsFactors = FALSE)
    net <- data.frame(source = c(rep("HIF1A",5), rep("EPAS1",4)),
                      target = c("a","b","c","d","e","a","b","c","f"),
                      stringsAsFactors = FALSE)
    inv <- data.frame(perturbed_gene = "HIF1A", n_up = 5000, n_down = 5000,
                      cell_line = "A", tissue = "T", dataset_id = "d1",
                      stringsAsFactors = FALSE)
    base <- mechanism_report(x, inventory = inv)
    full <- mechanism_report(x, inventory = inv, network = net)
    expect_equal(base$verdicts$class, full$verdicts$class)
    expect_equal(base$verdicts$verdict, full$verdicts$verdict)
    expect_equal(full$verdicts$class, "supported")
})

test_that("modules and warnings are printed once each, not once per pathway", {
    ## 同一件事按通路行重复：12 条通路会印 12 遍同一条警告；
    ## `{ATF6, XBP1}` 与 `{XBP1, ATF6}` 是同一组，不能印两次。
    x <- data.frame(pathway = paste0("p", 1:6), level = rep("B", 6),
                    verdict = rep("single-dataset", 6), reason = NA_character_,
                    driver = c("ATF6","XBP1","ATF6","XBP1","ATF6","XBP1"),
                    es = rep(-1, 6), padj = rep(0.01, 6), stringsAsFactors = FALSE)
    net <- data.frame(source = c(rep("ATF6", 6), rep("XBP1", 5)),
                      target = c("a","b","c","d","e","f","a","b","c","d","g"),
                      stringsAsFactors = FALSE)
    inv <- data.frame(perturbed_gene = c("ATF6","XBP1"),
                      n_up = c(2000, 600), n_down = c(2000, 500),
                      cell_line = c("A","B"), tissue = c("T","T"),
                      dataset_id = c("d1","d2"), stringsAsFactors = FALSE)
    out <- capture.output(print(mechanism_report(x, inventory = inv, network = net)))
    expect_equal(sum(grepl("^    \\{", out)), 1L)          # 只有一组
    expect_equal(sum(grepl("^    ATF6: available perturbations", out)), 1L)
})

test_that("support does not promise that future experiments are useless", {
    pe <- .make_pe("A", "perturbation-supported", NA, driver = "TF")
    inv <- data.frame(perturbed_gene = "TF", dataset_id = "d", cell_line = "A")
    out <- mechanism_report(pe, inventory = inv)
    expect_equal(out$verdicts$class, "supported")
    expect_false(any(grepl("No experiment would change|cannot be told apart|drags any pathway", out$summary)))
    expect_match(paste(out$summary, collapse = " "), "not identification")
    expect_match(out$next_steps$action, "independent perturbation")
})

test_that("all candidates retain supporting, conflicting and failed tests", {
    pe <- .make_pe("A", "perturbation-supported", NA, driver = "TF1")
    tests <- data.frame(pathway = "p", driver = c(rep("TF1", 4), "TF2", "TF3"),
        dataset = paste0("d", 1:6), cell_line = c("A", "B", "C", NA, "D", ""),
        es = c(-1, -1, 1, NA, 0.1, NA), padj = c(.01, .01, .01, NA, .8, NA),
        status = c("supports", "supports", "conflicts", "failed", "not-significant", "not-computable"),
        direction = "up", stringsAsFactors = FALSE)
    attr(pe, "perturbation_tests") <- tests
    inv <- data.frame(perturbed_gene = c("TF1", "TF2", "TF3"), dataset_id = paste0("i", 1:3), cell_line = "A")
    out <- mechanism_report(pe, inventory = inv)
    expect_identical(out$perturbation_tests, tests)
    expect_equal(out$candidate_evidence$driver, c("TF1", "TF2", "TF3"))
    expect_equal(out$candidate_evidence$n_support, c(2L, 0L, 0L))
    expect_equal(out$candidate_evidence$n_conflict, c(1L, 0L, 0L))
    expect_equal(out$candidate_evidence$n_failed, c(1L, 0L, 0L))
    expect_equal(out$candidate_evidence$assessment, c("mixed", "no-significant-support", "unavailable"))
    expect_equal(out$candidate_evidence$supporting_cells[1], "A, B")
    expect_match(paste(out$summary, collapse = " "), "both supporting and conflicting")
    expect_match(paste(capture.output(print(out)), collapse = " "), "support / conflict / failed")
    expect_setequal(out$gaps$driver, c("TF1", "TF2", "TF3"))
})

test_that("neighbourhood and breadth flags do not claim causal equivalence", {
    pe <- .make_pe("A", "perturbation-supported", NA, driver = "A")
    inv <- data.frame(perturbed_gene = c("A", "B"), dataset_id = c("d1", "d2"),
                      n_up = c(2000, 1), n_down = c(2000, 1))
    net <- data.frame(source = rep(c("A", "B"), each = 5), target = rep(letters[1:5], 2),
                      weight = rep(c(1, -1), each = 5))
    out <- mechanism_report(pe, drivers = c("A", "B"), inventory = inv, network = net)
    expect_equal(out$verdicts$module_size, 2L)
    expect_match(out$verdicts$global_effect_reason, "not proof of essentiality")
    expect_match(paste(out$summary, collapse = " "), "do not establish indistinguishability")
    expect_false(any(grepl("drags any pathway|cannot be told apart", out$summary)))
    expect_error(mechanism_report(pe, inventory = inv, network = net, min_overlap = NA_real_), "min_overlap")
    expect_error(mechanism_report(pe, inventory = inv, global_effect_de = -1), "global_effect_de")
})

test_that("breadth lookup resolves driver aliases", {
    inv <- data.frame(perturbed_gene = "NR3C1", n_up = 2000, n_down = 2000)
    x <- clusterProfiler:::.mechanism_global_effect(inv, c("GR", "NR3C1"))
    expect_equal(unname(x), c(4000, 4000))
})

test_that("unknown inventories cannot be reported as globally absent", {
    inv <- data.frame(perturbed_gene = character(), dataset_id = character(), cell_line = character())
    attr(inv, "query_status") <- data.frame(driver = c("A", "B"), status = c("failed", "empty"), message = c("outage", ""))
    pe <- .make_pe("C", "unsupported", "no-dataset")
    out <- mechanism_report(pe, drivers = c("A", "B"), inventory = inv)
    expect_true(is.na(out$gaps$n_datasets[1]))
    expect_equal(out$gaps$n_datasets[2], 0L)
    expect_match(out$gaps$gap[1], "coverage unknown")
    expect_match(out$gaps$gap[2], "resolved by this query")
    attr(inv, "query_status") <- NULL
    supplied <- mechanism_report(pe, drivers = "A", inventory = inv[0, ])
    expect_match(supplied$gaps$gap, "supplied inventory")
})

test_that("unresolved requested drivers remain visible in report gaps", {
    pe <- .make_pe("C", "unsupported", "no-dataset")
    attr(pe, "drivers") <- c("TF1", "TF2")
    attr(pe, "perturbation_tests") <- clusterProfiler:::.pathway_perturbation_tests_empty()
    inv <- data.frame(perturbed_gene = character(), dataset_id = character(), cell_line = character())
    out <- mechanism_report(pe, inventory = inv)
    expect_equal(out$gaps$driver, c("TF1", "TF2"))
    expect_equal(out$gaps$n_datasets, c(0L, 0L))
    expect_true(all(grepl("supplied inventory", out$gaps$gap)))
})

test_that("a row subset cannot import unrelated pathway conflicts", {
    pe <- .make_pe("A", "perturbation-supported", NA, driver = "TF1", pathway = "p1")
    attr(pe, "perturbation_tests") <- data.frame(pathway = c("p1", "p2", "p2"),
        driver = c("TF1", "TF2", "TF2"), dataset = c("d1", "d2", "d3"),
        cell_line = c("A", "B", "C"), status = c("supports", "supports", "conflicts"))
    inv <- data.frame(perturbed_gene = c("TF1", "TF2"))
    out <- mechanism_report(pe, inventory = inv)
    expect_equal(out$candidate_evidence$driver, "TF1")
    expect_equal(out$perturbation_tests$pathway, "p1")
    expect_false(any(grepl("both supporting and conflicting", out$summary)))
})

test_that("empty report exposes stable typed evidence schemas", {
    out <- mechanism_report(data.frame())
    expect_equal(nrow(out$candidate_evidence), 0L)
    expect_true(all(c("pathway", "driver", "n_support", "n_conflict", "assessment",
        "n_noncomputable", "n_nonsignificant", "n_attempted") %in% names(out$candidate_evidence)))
    expect_type(out$verdicts$evidence_status, "character")
    expect_type(out$verdicts$n_computable, "integer")
    expect_type(out$verdicts$incomplete, "logical")
    expect_false(any(grepl("No experiment would change", out$summary)))
})

test_that("ledger states lead the report without reclassifying legacy verdicts", {
    paths <- c("allNC", "NSonly", "conflictB", "supportB", "mixedA", "partialfailure")
    pe <- do.call(rbind, lapply(seq_along(paths), function(i)
        .make_pe(if (i == 5L) "A" else "B", "single-dataset", "low-support",
                 driver = "TF", pathway = paths[i])))
    ledger <- data.frame(pathway = rep(paths, c(2, 1, 1, 1, 2, 3)), driver = "TF",
        dataset = paste0("d", 1:10), cell_line = "A",
        status = c("not-computable", "not-computable", "not-significant", "conflicts",
            "supports", "supports", "conflicts", "supports", "not-significant", "failed"))
    attr(pe, "perturbation_tests") <- ledger
    out <- mechanism_report(pe, inventory = data.frame(perturbed_gene = "TF"))
    expect_equal(out$verdicts$evidence_status,
        c("noncomputable", "undetermined", "conflict", "support", "mixed", "support"))
    expect_equal(out$verdicts$class, c(rep("undetermined", 4), "supported", "undetermined"))
    expect_identical(out$verdicts$verdict, pe$verdict)
    expect_identical(out$verdicts$level, pe$level)
    expect_identical(out$perturbation_tests, ledger)
    expect_equal(out$verdicts$n_computable, c(0L, 1L, 1L, 1L, 2L, 2L))
    expect_equal(out$verdicts$n_attempted, c(2L, 1L, 1L, 1L, 2L, 3L))
    expect_equal(out$verdicts$availability, c(rep("available", 5), "partial"))
    expect_equal(out$verdicts$incomplete, c(TRUE, rep(FALSE, 4), TRUE))
    expect_equal(out$candidate_evidence$n_noncomputable, c(2L, rep(0L, 5)))
    expect_equal(out$candidate_evidence$n_nonsignificant, c(0L, 1L, 0L, 0L, 0L, 1L))
    expect_equal(out$candidate_evidence$assessment,
        c("unavailable", "no-significant-support", "conflict-only", "support-only", "mixed", "support-only"))
    printed <- capture.output(print(out, n = 10))
    expect_equal(printed[grepl("^  (support|conflict|mixed|undetermined|noncomputable):", printed)],
        paste0("  ", c("support", "conflict", "mixed", "undetermined", "noncomputable"), ":"))
    expect_match(paste(printed, collapse = "\n"), "2/3 computable, 0 noncomputable, 1 failed \\[partial\\]")
    expect_match(paste(out$summary, collapse = " "), "not.*grading pass")
})

test_that("missing ledgers and lookup outages cannot invent support or noncomputability", {
    pe <- rbind(.make_pe("A", "perturbation-supported", NA, pathway = "a"),
                .make_pe("B", "direction-conflict", "direction-conflict", pathway = "b"))
    inv <- data.frame(perturbed_gene = character())
    out <- mechanism_report(pe, inventory = inv)
    expect_equal(out$verdicts$evidence_status, rep("undetermined", 2))
    expect_equal(out$verdicts$availability, rep("missing-ledger", 2))
    expect_equal(out$verdicts$class, c("supported", "conflicting"))
    expect_true(all(out$verdicts$incomplete))
    attr(pe, "perturbation_tests") <- clusterProfiler:::.pathway_perturbation_tests_empty()
    for (status in c("failed", "unavailable", "not-run", "empty")) {
        attr(pe, "dataset_lookup") <- data.frame(driver = "TF", status = status)
        out <- mechanism_report(pe, inventory = inv)
        expect_equal(out$verdicts$evidence_status, rep("undetermined", 2))
        expect_equal(out$verdicts$availability, rep(switch(status,
            failed = "lookup-unavailable", unavailable = "lookup-unavailable",
            `not-run` = "unqueried", empty = "no-dataset"), 2))
        expect_equal(out$verdicts$n_attempted, c(0L, 0L))
    }
    attr(pe, "perturbation_tests") <- data.frame(pathway = "a", driver = "TF",
        dataset = "d", cell_line = "A", status = "failed")
    attr(pe, "dataset_lookup") <- NULL
    out <- mechanism_report(pe, inventory = inv)
    expect_equal(out$verdicts$evidence_status, c("undetermined", "undetermined"))
    expect_equal(out$verdicts$availability, c("failed", "not-tested"))
    expect_equal(out$verdicts$n_failed, c(1L, 0L))
    attr(pe, "perturbation_tests") <- data.frame(pathway = "a", driver = "TF",
        dataset = "d", cell_line = "A", status = "supports")
    for (status in c("empty", "failed", "not-run")) {
        attr(pe, "dataset_lookup") <- data.frame(driver = c("TF", "missing"),
            status = c("available", status))
        out <- mechanism_report(pe, inventory = inv)
        expect_equal(out$verdicts$evidence_status[1], "support")
        expect_equal(out$verdicts$availability[1], "partial")
        expect_true(out$verdicts$incomplete[1])
        expect_equal(out$verdicts$n_failed[1], 0L)
    }
    attr(pe, "dataset_lookup") <- NULL
    attr(pe, "perturbation_tests")$status <- NA_character_
    unknown <- clusterProfiler:::.pathway_evidence_status(pe)
    expect_equal(unknown$availability[1], "unknown")
    expect_equal(unknown$evidence_status[1], "undetermined")
})

test_that("unknown and missing statuses never enter computable candidate denominators", {
    pe <- .make_pe("B", "unsupported", NA, driver = "TF")
    attr(pe, "perturbation_tests") <- data.frame(pathway = "p", driver = "TF",
        dataset = paste0("d", 1:5), cell_line = "A",
        status = c("supports", "not-significant", NA, "unexpected", "not-computable"))
    out <- mechanism_report(pe, inventory = data.frame(perturbed_gene = "TF"))
    expect_equal(out$candidate_evidence$n_tests, 4L) # Historical definition is preserved.
    expect_equal(out$candidate_evidence$n_computable, 2L)
    expect_identical(out$candidate_evidence$n_computable, out$verdicts$n_computable)
    expect_equal(out$candidate_evidence$n_attempted, 5L)
    printed <- capture.output(print(out))
    expect_match(printed[grepl("p <- TF:", printed)], "2/5 computable")
    expect_false(any(grepl("4/5 computable", printed)))
    empty <- clusterProfiler:::.mechanism_candidate_evidence(NULL)
    expect_type(empty$n_computable, "integer")
})

test_that("serialized legacy reports keep pathways and candidate rows visible", {
    pe <- .make_pe("A", "perturbation-supported", NA, driver = "TF", pathway = "saved-pathway")
    attr(pe, "perturbation_tests") <- data.frame(pathway = "saved-pathway", driver = "TF",
        dataset = "d", cell_line = "A", status = "supports")
    old <- mechanism_report(pe, inventory = data.frame(perturbed_gene = "TF"))
    old$verdicts <- old$verdicts[, c("pathway", "class", "verdict", "reason", "detail",
        "side", "action", "module", "module_size", "global_effect_reason")]
    old$candidate_evidence <- old$candidate_evidence[, c("pathway", "driver", "n_tests",
        "n_support", "n_conflict", "n_failed", "supporting_cells", "conflicting_cells", "assessment")]
    for (has_ledger in c(TRUE, FALSE)) {
        if (!has_ledger) old$perturbation_tests <- NULL
        saved <- unserialize(serialize(old, NULL))
        before <- serialize(saved, NULL)
        printed <- capture.output(returned <- print(saved))
        expect_match(paste(printed, collapse = "\n"), "legacy saved verdicts")
        expect_match(paste(printed, collapse = "\n"), "evidence status and grade not inferred")
        expect_true(any(grepl("- saved-pathway", printed, fixed = TRUE)))
        expect_true(any(grepl("saved-pathway <- TF: 1 / 0 / 0", printed, fixed = TRUE)))
        expect_true(any(grepl("computable denominator unavailable", printed, fixed = TRUE)))
        expect_false(any(grepl("^  (support|conflict|mixed|undetermined|noncomputable):", printed)))
        expect_identical(serialize(saved, NULL), before)
        expect_identical(returned, saved)
    }
})
