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
                 "no-direction", "no-fold-changes", "no-test-completed",
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

test_that("verdict classes split supported, excluded and undetermined", {
    pe <- rbind(.make_pe("A", "perturbation-supported", NA),
                .make_pe("B", "direction-conflict", "direction-conflict"),
                .make_pe("C", "unsupported", "too-small"))
    expect_equal(clusterProfiler:::.mechanism_verdict_class(pe),
                 c("supported", "excluded", "undetermined"))
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
    expect_match(g$gap[1], "no dataset from Airway")
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
    expect_equal(sum(rep$verdicts$class == "undetermined"), 2L)
    expect_equal(sum(rep$verdicts$class == "excluded"), 1L)
})

test_that("the report says plainly when nothing would change the verdicts", {
    pe <- .make_pe("C", "coherent-refuted", "not-significant")
    rep <- mechanism_report(pe, drivers = character(0))
    expect_equal(nrow(rep$next_steps), 0L)
    expect_match(rep$summary[length(rep$summary)], "No experiment")
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
    expect_match(paste(rep$summary, collapse = " "), "No experiment would change")
    rep2 <- mechanism_report(pe, drivers = character(0))
    expect_gt(nrow(rep2$next_steps), 0L)
})

test_that("print covers the three sections and the inventory", {
    pe <- rbind(.make_pe("C", "unsupported", "single-cell-line", driver = "GR", pathway = "a"),
                .make_pe("C", "coherent-refuted", "not-significant", pathway = "b"))
    rep <- mechanism_report(pe, inventory = data.frame(
        perturbed_gene = "GR", dataset_id = "D1", cell_line = "HEK293T",
        tissue = "Kidney", n_samples = 3, stringsAsFactors = FALSE))
    out <- capture.output(print(rep))
    expect_true(any(grepl("mechanism report", out)))
    expect_true(any(grepl("excluded", out)))
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
    ## 旁系同源 / 同复合物 / 并行支路本来就不能用表达实验分开。
    ## 报一个模块比给一个假的精确排名诚实。
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
    expect_equal(sum(grepl("^    ATF6: knocking down", out)), 1L)
})
