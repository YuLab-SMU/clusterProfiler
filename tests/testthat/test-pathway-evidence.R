## 通路模式的第一批测试。
##
## 这些测的是**判断逻辑**，不是网络。`pathway_evidence()` 在有驱动者时会调 GPSAdb，
## 所以涉及网络的用例一律 skip；离线用例覆盖方向推导、成员一致性、分级与陈述。

## 把若干行拼成一个富集结果，用于"一条匹配、一条不匹配"的混合输入。
rbind_result <- function(...) {
    parts <- lapply(list(...), function(z) as.data.frame(z, stringsAsFactors = FALSE))
    do.call(rbind, parts)
}

pathway_fixture_result <- function() {
    data.frame(
        ID = c("GO:1", "GO:2", "GO:3"),
        Description = c("response to hypoxia", "coherent pathway", "incoherent pathway"),
        geneID = c("A/B/C/D", "E/F/G/H", "I/J/K/L"),
        Count = 4L, p.adjust = c(1e-10, 1e-5, 1e-3),
        stringsAsFactors = FALSE)
}

test_that("direction is derived from the members when the result has none", {
    ## ORA 结果本身不带方向：它只说通路富集了，不说往哪边动。
    ## 第一版这里返回 NULL，导致每一行都变成 untestable——这是跑起来才发现的。
    r <- pathway_fixture_result()[1, ]
    members <- c("A", "B", "C", "D")
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5)
    expect_identical(.pathway_direction(r, members, signed), "up")
    expect_identical(.pathway_direction(r, members, -signed), "down")
    ## 没有 signed 时不猜。
    expect_null(.pathway_direction(r, members, NULL))
})

test_that("member support measures the fraction moving together", {
    ## 少数基因驱动的高分不能算整条通路的支持。
    members <- c("A", "B", "C", "D")
    expect_equal(.pathway_support(members, c(A = 1, B = 2, C = 3, D = 4), "up"), 1)
    expect_equal(.pathway_support(members, c(A = 1, B = -2, C = 3, D = -4), "up"), 0.5)
    ## 成员一个都不在被检验的基因里时，不是 0，是 NA——"没测"不等于"没通过"。
    expect_true(is.na(.pathway_support(members, c(Z = 1), "up")))
})

test_that("no signed input is untestable, and says which", {
    ## 没有 signed 是**没测**，不是测了不成立。
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, max_pathways = 3)
    expect_equal(nrow(out), 3)
    expect_true(all(out$level == "C"))
    expect_true(all(out$verdict == "untestable"))
    expect_true(all(grepl("No signed fold changes were supplied", out$statement)))
})

test_that("support derived from the members themselves cannot promote to level B", {
    ## 审核发现的缺陷：方向由成员多数推导时，support 恒 >= 0.5，`min_support = 0.5`
    ## 永远成立，于是每条有方向的通路都白得一个 B。方向自己来自成员，再拿成员的
    ## 一致性当证据，是同一批基因既定方向又证方向。
    r <- pathway_fixture_result()[1:2, ]
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = -2, F = -1, G = -3, H = -0.5)
    out <- pathway_evidence(r, signed = signed, max_pathways = 2)
    expect_true(all(out$level == "C"))
    expect_false(any(out$verdict == "coherent-untested"))
    ## 但 support 仍然照实报出来，它不是错，只是不能当独立证据。
    expect_true(all(out$support == 1))
})

test_that("an externally supplied direction makes member support meaningful again", {
    ## 方向来自结果本身（外部信息）时，成员一致性才是独立证据，可以给 B。
    r <- pathway_fixture_result()[1, ]
    r$pathway_direction <- "up"
    signed <- c(A = 2, B = 1.5, C = 3, D = -0.5)
    out <- pathway_evidence(r, signed = signed)
    expect_equal(out$level, "B")
    expect_equal(out$verdict, "coherent-untested")
    expect_equal(out$support, 0.75)
})

test_that("members splitting evenly are unsupported, not untested", {
    ## 给了 signed、成员却对半分：这是**测了，没有共同方向**，与"没测"是两个结论。
    r <- pathway_fixture_result()[3, ]
    signed <- c(I = 2, J = -1, K = 3, L = -4)
    out <- pathway_evidence(r, signed = signed)
    expect_equal(out$level, "C")
    expect_equal(out$verdict, "unsupported")
    expect_match(out$statement, "split evenly")
})

test_that("an empty result returns an empty frame with the right columns", {
    out <- pathway_evidence(pathway_fixture_result()[0, ])
    expect_equal(nrow(out), 0)
    expect_true(all(c("pathway", "level", "statement") %in% names(out)))
})

test_that("level A requires two cell lines, and one dataset is not enough", {
    ## 实测依据：四个驱动者里三个换一个细胞系就不成立，效应量符号也翻转
    ## （results/stage3_refutation.tsv）。只留最强一份的实现让 A 级等于
    ## "某一份数据支持"，这个含义太弱，所以 A 级现在要求复现。
    ##
    ## 这里只能测判断逻辑，不调网络：用 skip 保护，真值由 pathway_evidence 内部算。
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))
    r <- pathway_fixture_result()[1, ]
    out <- pathway_evidence(r, drivers = "HIF1A", signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_true("n_cell_lines_supported" %in% names(out))
    expect_true("n_tests_run" %in% names(out))
    ## 若有显著支持但只有单一细胞系，不能是 A。
    if (out$n_tests_run > 0 && out$padj < 0.05 && out$n_cell_lines_supported < 2) {
        expect_false(out$level == "A")
        expect_equal(out$verdict, "single-dataset")
    }
})

test_that("every row says why it stopped short of level A", {
    ## `reason` 存在的理由：`single-dataset` 与 `direction-conflict` 都落在 B，
    ## 一个是证据不够多、一个是证据指向相反；三种 C（成员找不到、通路太小、
    ## 测了不显著）也需要分开。否则读者无法区分"没找到"与"没法找"。
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, drivers = NULL, signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_true("reason" %in% names(out))
    ## 没有驱动者时，B 级只能来自"通路动了但没有扰动检验"。
    if (any(out$level == "B")) {
        expect_true(all(out$reason[out$level == "B"] %in% "no-perturbation-test"))
    }
    ## A 级不该带降级原因。
    expect_true(all(is.na(out$reason[out$level == "A"])))
})

test_that("the result carries a class so print.pathway_evidence is used", {
    ## 实测这份代码写了之后一次都没被派发过：没有 class，返回裸 data.frame，
    ## 用户看到 20 列原始表格。是最难发现的一类"看起来有、其实没跑"。
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, drivers = NULL, signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_s3_class(out, "pathway_evidence")
    expect_true(is.data.frame(out))
})

test_that("the statement agrees with the numbers printed beside it", {
    ## 审核发现的缺陷，而且是**修过一次仍然存在**的那种：`why_c` 把
    ## "成员方向分散"与"方向由成员自己推导"合成一个 `low-support`，
    ## 于是 support = 100% 的通路印出 "100% of members moved up, below the
    ## support floor of 50%"。上一轮加的那条 `else` 分支永远走不到。
    ##
    ## 这条测试的价值在于它**断言两列之间的关系**，而不只是各自的取值——
    ## 只测 level 和 support 的话，这个矛盾照样活下来。
    ## 需要 >= 5 个成员，否则先命中 `too-small`，测不到这一支。
    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)
    out <- pathway_evidence(r, signed = signed)

    expect_equal(out$support, 1)
    expect_false(grepl("below the support floor", out$statement))
    expect_match(out$statement, "inferred from these same members")
    expect_equal(out$reason, "member-derived")

    ## 反向：支撑率真的低于下限时，才可以说"低于下限"。
    r2 <- r
    r2$pathway_direction <- "up"          # 方向来自外部，support 才是独立证据
    signed2 <- c(A = 2, B = -1.5, C = -3, D = -0.5, E = -1)
    out2 <- pathway_evidence(r2, signed = signed2)
    expect_lt(out2$support, 0.5)
    expect_match(out2$statement, "below the support floor")
    expect_equal(out2$reason, "low-support")
})

test_that("max_pathways = 0 does not crash", {
    ## 原来 `head(d, 0)` 之后没有守卫，do.call(rbind, list()) 返回 NULL，
    ## 再 attr(out, ...) <- 就报 "attempt to set an attribute on NULL"。
    r <- pathway_fixture_result()
    out <- expect_no_error(pathway_evidence(r, max_pathways = 0L,
                                            signed = c(A = 2, B = 1.5)))
    ## 夹到 1 而不是报错，也不是返回 0 行——问"给我零条"没有有意义的答案。
    expect_equal(nrow(out), 1)
    expect_s3_class(out, "pathway_evidence")
})

test_that("an empty result carries the same class as a full one", {
    ## 两种返回形态不一致时，空结果不走 print 方法，用户在同一函数上
    ## 看到两种输出。
    expect_s3_class(pathway_evidence(pathway_fixture_result()[0, ]), "pathway_evidence")
})

test_that("print shows a placeholder rather than NA for untestable rows", {
    ## `%||%` 只接 NULL，接不住 NA，原来方向印成 "NA"、支撑率印成 "NA%"。
    r <- pathway_fixture_result()[1, ]
    out <- pathway_evidence(r)          # 没有 signed -> 方向是 NA
    txt <- paste(utils::capture.output(print(out)), collapse = "\n")
    expect_false(grepl("NA%", txt, fixed = TRUE))
    expect_match(txt, "\\?")
})

test_that("a failed perturbation call is reported as an outage, not as weak evidence", {
    ## 审核发现：`ask()` 把调用失败和"调用成功但没结果"都当成 next，
    ## 于是一次断网让 `tested` 变空，再被解释成"成员方向不一致"。
    ## 实测 HIF1A（4 个细胞系）在模拟断网下拿到 `low-support`、零报警。
    ##
    ## 需要网络先把驱动者解析成数据集，所以走同一道网络开关。
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)

    ns <- asNamespace("gpsadb")
    orig <- get("gpsadb_test_gene_set", envir = ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", orig, ns = "gpsadb"), add = TRUE)
    assignInNamespace("gpsadb_test_gene_set",
                      function(...) stop("simulated outage"), ns = "gpsadb")

    expect_warning(out <- pathway_evidence(r, drivers = "HIF1A", signed = signed),
                   "GPSAdb calls failed")
    expect_equal(out$reason, "no-test-completed")
    expect_equal(out$n_tests_run, 0)
    expect_match(out$statement, "outage")
})

test_that("subsetting columns drops the class so print keeps working", {
    ## `[.data.frame` 会把自定义 class 一起带过去，于是 `out[, c("pathway","level")]`
    ## 仍声称自己是 pathway_evidence，而 print 要读 driver/reason/support 这些
    ## 已经不存在的列，`is.na(NULL)` 是 logical(0)，`if` 收到长度 0 就报
    ## "argument is of length zero"。按列取子集太常见了，必须处理。
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, signed = c(A = 2, B = 1.5, C = 3, D = 0.5))

    sub <- out[, c("pathway", "level")]
    expect_false(inherits(sub, "pathway_evidence"))
    expect_no_error(print(sub))

    ## 行子集保留全部列，class 应当留着。
    expect_s3_class(out[1, ], "pathway_evidence")
    ## 单列取子集是向量，不该被强加 class。
    expect_false(inherits(out[, "level"], "pathway_evidence"))
})

test_that("empty separators and NA are not counted as members", {
    ## `"A/B/C/D/E/"` 原来数出 6 个成员，`n_members` 与陈述里的数字都多 1；
    ## `geneID = NA` 会得到一个叫 "NA" 的成员，把"没有成员"报成 1 个。
    base <- pathway_fixture_result()[1, ]
    mk <- function(g) { base$geneID <- g; base }

    members <- function(g) clusterProfiler:::.pathway_members(mk(g))
    expect_length(members("A/B/C/D/E"), 5)
    expect_length(members(" A/B/C/D/E/ "), 5)
    expect_length(members("A//B/C/D/E"), 5)
    expect_length(clusterProfiler:::.pathway_members(mk(NA_character_)), 0)
})

test_that("no common direction names which of the three cases it is", {
    ## 没给 signed / 成员不在差异表里 / 成员真的对半分 —— 三个不同结论。
    ## 原来后两者共用 "Members split evenly"，一条 0 成员的通路也说成"对半分"。
    base <- pathway_fixture_result()[1, ]

    r1 <- base
    expect_equal(pathway_evidence(r1)$reason, "no-fold-changes")

    r2 <- base
    r2$geneID <- "ZZZ/YYY"
    r2$Description <- "members not in the differential table"
    out2 <- pathway_evidence(rbind_result(r2, base), signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_equal(out2$reason[out2$pathway == r2$Description], "no-members")
    expect_match(out2$statement[out2$pathway == r2$Description],
                 "appear in the signed gene list")

    r3 <- base
    r3$geneID <- "A/B"
    r3$pathway_direction <- NA
    out3 <- pathway_evidence(r3, signed = c(A = 2, B = -1))
    expect_equal(out3$reason, "no-direction")
    expect_match(out3$statement, "split evenly")
})


test_that("the documented evidence() pathway table keeps its names", {
    ## 文档允许 `x` 是 `evidence()` bundle 的 `$pathways` 表，而那张表用的是
    ## `pathway` / `pathway_description` 两个列名，不是 `ID` / `Description`。
    ## 只认后者时每行通路名都是 NA——不报错，只是名字空了。
    pw <- data.frame(
        pathway = "R-HSA-1234176",
        pathway_description = "Oxygen-dependent proline hydroxylation of HIF Alpha",
        pathway_direction = "up", n_pathway_genes = 5L, pathway_padj = 1e-8,
        regulator = "VHL", source = "reactome", term = "R-HSA-1234176",
        Count = 5L, setSize = 40L, padj_within_pathway = 1e-5,
        geneID = "A/B/C/D/E", stringsAsFactors = FALSE)

    out <- pathway_evidence(pw, signed = c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1))
    expect_false(is.na(out$pathway))
    expect_equal(out$pathway, "Oxygen-dependent proline hydroxylation of HIF Alpha")
})

test_that("a GSEA result takes its direction from NES, which beats the members", {
    ## `gseaResult` 自带 NES，方向应当来自它；成员多数方向只在没有 NES 时才用。
    g <- data.frame(ID = "GO:1", Description = "response to hypoxia",
                    setSize = 5L, enrichmentScore = 0.8, NES = -2.0,
                    pvalue = 0.001, p.adjust = 0.01, qvalue = 0.01,
                    rank = 100, leading_edge = "tags=50%",
                    core_enrichment = "A/B/C/D/E", geneID = "A/B/C/D/E",
                    stringsAsFactors = FALSE)
    ## 成员全部上调，但 NES 说这条通路是 down —— 方向应当听 NES。
    out <- pathway_evidence(g, signed = c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1))
    expect_equal(out$direction, "down")
    ## 方向与成员相反时，支持率是 0，不能当证据。
    expect_equal(out$support, 0)
    expect_equal(out$reason, "low-support")
})

test_that("an outage is not hidden behind 'no perturbation test identified a driver'", {
    ## 方向来自外部且支撑率够时，B 分支会先命中，把"检验全失败"说成
    ## "没有检验能指出驱动者"。实测断网就是这句话，看不出是网络问题。
    ##
    ## 需要网络解析数据集，走同一道开关。
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)

    ns <- asNamespace("gpsadb")
    orig <- get("gpsadb_test_gene_set", envir = ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", orig, ns = "gpsadb"), add = TRUE)
    assignInNamespace("gpsadb_test_gene_set",
                      function(...) stop("simulated outage"), ns = "gpsadb")

    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A", signed = signed))
    expect_equal(out$reason, "no-test-completed")
    expect_equal(out$level, "C")
})

test_that("a driver that was tested and refuted is not called untested", {
    ## 检验跑了、全不显著，与根本没跑是两件事。原来两者共用
    ## `coherent-untested`，一个跑了 4 个数据集全不显著的通路会说成"没测"——
    ## `untested` 这个词本身就是错的。
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)

    ns <- asNamespace("gpsadb")
    orig <- get("gpsadb_test_gene_set", envir = ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", orig, ns = "gpsadb"), add = TRUE)
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 0.5, es = -0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")

    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A", signed = signed))
    expect_gt(out$n_tests_run, 0)
    expect_equal(out$verdict, "coherent-refuted")
    expect_equal(out$reason, "not-significant")
    expect_match(out$statement, "found no significant movement")

    ## 没给驱动者时仍然是"没测"，两者不能混。
    out2 <- pathway_evidence(r, signed = signed)
    expect_equal(out2$verdict, "coherent-untested")
    expect_equal(out2$reason, "no-perturbation-test")
})

test_that("the level-A statement admits what it does not establish", {
    ## 实测：同一条缺氧通路上，从库里随机取 30 个驱动者，有 4 个（全是必需基因）
    ## 同样拿到 A。所以 A 只说明"敲掉它这条路会动"，不说明"它就是驱动者"。
    ## 陈述里必须有这句话，否则读者会把"一致"读成"认定"。
    ##
    ## 打桩造出一个真正的 A 级（两份细胞系、同向、显著），不依赖网络。
    ## 先前写的版本用 3 个成员的通路，永远拿不到 A，于是整条测试被 skip 掉了——
    ## 跳过等于没有覆盖，所以改成打桩。
    ns <- asNamespace("gpsadb")
    orig_s <- get("gpsadb_search_datasets", envir = ns)
    orig_t <- get("gpsadb_test_gene_set", envir = ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", orig_s, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", orig_t, ns = "gpsadb")
    }, add = TRUE)

    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        list(data.frame(perturbed_gene = rep(perturbed_gene, 2),
                        method = "knockdown",
                        cell_line = c("HeLa", "A-498"),
                        n_samples = c(4L, 4L), in_reference = c(TRUE, TRUE),
                        dataset_id = c("D001", "D002"), stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 1e-6, es = -0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"          # 方向来自外部，es < 0 即 activator
    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A",
                                             signed = c(A = 2, B = 1.5, C = 3,
                                                        D = 0.5, E = 1)))

    expect_equal(out$level, "A")
    expect_match(out$statement, "Consistent with drive")
    expect_match(out$statement, "not evidence that")
})

test_that("a single-cell-line row says the cap is the data, not the evidence", {
    ## B 级里有两类：方向冲突（真的相反）与单细胞系（只是不够）。
    ## 后者实测在 11 个已知真值的数据集上指向全对，措辞不能把它写成次品。
    ns <- asNamespace("gpsadb")
    orig_s <- get("gpsadb_search_datasets", envir = ns)
    orig_t <- get("gpsadb_test_gene_set", envir = ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", orig_s, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", orig_t, ns = "gpsadb")
    }, add = TRUE)
    ## 只给一个细胞系——这正是"够不到 A"的那种数据。
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        list(data.frame(perturbed_gene = perturbed_gene, method = "knockdown",
                        cell_line = "HeLa", n_samples = 4L, in_reference = TRUE,
                        dataset_id = "D001", stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 1e-6, es = -0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A",
                                             signed = c(A = 2, B = 1.5, C = 3, D = 2, E = 1)))
    expect_equal(out$level[1], "B")
    expect_equal(out$reason[1], "single-cell-line")
    expect_match(out$statement[1], "limit of the available data")
    expect_match(out$statement[1], "consistent perturbation evidence, but under-evidenced")
})

test_that("a direction conflict is not described as thin data", {
    ## 两类 B 的措辞必须能分开：方向冲突是证据反对，单细胞系是数据不够。
    ns <- asNamespace("gpsadb")
    orig_s <- get("gpsadb_search_datasets", envir = ns)
    orig_t <- get("gpsadb_test_gene_set", envir = ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", orig_s, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", orig_t, ns = "gpsadb")
    }, add = TRUE)
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        list(data.frame(perturbed_gene = perturbed_gene, method = "knockdown",
                        cell_line = "HeLa", n_samples = 4L, in_reference = TRUE,
                        dataset_id = "D001", stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    ## 通路 up、es 也 > 0 —— 敲低把它推得更高，读作压制而不是驱动。
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 1e-6, es = 0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A",
                                             signed = c(A = 2, B = 1.5, C = 3, D = 2, E = 1)))
    expect_equal(out$verdict[1], "direction-conflict")
    expect_match(out$statement[1], "repression")
    expect_false(grepl("limit of the available data", out$statement[1]))
})

## Namespace stubs are restored even when an assertion fails; no network calls.
with_pathway_stubs <- function(datasets, outcomes, code) {
    ns <- asNamespace("gpsadb")
    old_search <- get("gpsadb_search_datasets", ns)
    old_test <- get("gpsadb_test_gene_set", ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", old_search, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", old_test, ns = "gpsadb")
    }, add = TRUE)
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        z <- datasets[datasets$driver == perturbed_gene, , drop = FALSE]
        z$n_samples <- rep(4L, nrow(z))
        z$in_reference <- rep(TRUE, nrow(z))
        list(z)
    }, ns = "gpsadb")
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        z <- outcomes[[id]]
        if (inherits(z, "error")) stop(z)
        list(z)
    }, ns = "gpsadb")
    force(code)
}

pathway_test_input <- function() {
    data.frame(Description = "programme", geneID = "A/B/C/D/E",
               pathway_direction = "up", stringsAsFactors = FALSE)
}

pathway_test_signed <- c(A = 1, B = 2, C = 3, D = 4, E = 5)

test_that("replication never pools cell lines across candidates", {
    ds <- data.frame(driver = c("TF1", "TF2"), dataset_id = c("D1", "D2"),
                     cell_line = c("cellA", "cellB"))
    got <- list(D1 = data.frame(es = -1, padj = 0.001),
                D2 = data.frame(es = -1, padj = 0.002))
    with_pathway_stubs(ds, got, {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(),
            drivers = c("TF1", "TF2"), signed = pathway_test_signed))
        expect_equal(out$level, "B")
        expect_equal(out$n_cell_lines_supported, 1L)
        expect_equal(out$driver, "TF1")
        ledger <- attr(out, "perturbation_tests")
        expect_equal(ledger$driver, c("TF1", "TF2"))
        expect_equal(ledger$status, c("supports", "supports"))
    })
})

test_that("significance cannot bypass measurable member support", {
    ds <- data.frame(driver = "TF1", dataset_id = c("D1", "D2"),
                     cell_line = c("cellA", "cellB"))
    got <- list(D1 = data.frame(es = -1, padj = 0.001),
                D2 = data.frame(es = -1, padj = 0.002))
    with_pathway_stubs(ds, got, {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = c(A = 1, B = -2, C = -3, D = -4, E = -5)))
        expect_equal(out$level, "C")
        expect_equal(out$reason, "low-support")
        expect_match(out$statement, "below the support floor")
        expect_equal(attr(out, "perturbation_tests")$status, rep("supports", 2))
        inferred <- pathway_test_input()
        inferred$pathway_direction <- NULL
        strict <- suppressWarnings(pathway_evidence(inferred, drivers = "TF1",
            min_support = 0.8, signed = c(A = 1, B = 2, C = 3, D = -4, E = -5)))
        expect_equal(strict$level, "C")
        expect_equal(strict$reason, "low-support")
        unknown <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1"))
        expect_equal(unknown$level, "C")
        ## Direction is given, but no fold changes: the missing input is the
        ## measurements, not the identifiers.
        expect_equal(unknown$reason, "no-fold-changes")
    })
})

test_that("the ledger retains all outcomes and warns about partial failures", {
    ds <- data.frame(driver = "TF1", dataset_id = paste0("D", 1:6),
                     cell_line = paste0("cell", 1:6))
    got <- list(D1 = data.frame(es = -1, padj = 0.001),
                D2 = data.frame(es = -1, padj = 0.002),
                D3 = data.frame(es = 1, padj = 0.003),
                D4 = data.frame(es = -1, padj = 0.5),
                D5 = simpleError("offline failure"),
                D6 = data.frame(es = NA_real_, padj = NA_real_))
    with_pathway_stubs(ds, got, {
        expect_warning(out <- pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = pathway_test_signed), "1 of 6 GPSAdb calls failed")
        ledger <- attr(out, "perturbation_tests")
        expect_identical(ledger$status, c("supports", "supports", "conflicts",
                                         "not-significant", "failed", "not-computable"))
        expect_identical(ledger$direction, rep("up", 6))
        expect_type(ledger$es, "double")
        expect_type(ledger$padj, "double")
        expect_equal(out$n_tests_run, 4L)
        expect_equal(out$n_cell_lines_supported, 2L)
        expect_equal(out$level, "A")
        expect_match(out$statement, "contradictory across datasets")
        expect_equal(attr(out, "network_calls"), 6L)
    })
})

test_that("empty and untestable results retain a typed ledger", {
    empty <- suppressWarnings(pathway_evidence(pathway_test_input()[0, ]))
    ledger <- attr(empty, "perturbation_tests")
    expect_identical(names(ledger), c("pathway", "driver", "dataset", "cell_line",
                                     "es", "padj", "status", "direction"))
    expect_equal(nrow(ledger), 0L)
    expect_type(ledger$es, "double")
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = simpleError("must not call")), {
        x <- pathway_test_input()
        x$pathway_direction <- NULL
        out <- suppressWarnings(pathway_evidence(x, drivers = "TF1"))
        expect_equal(attr(out, "network_calls"), 0L)
        expect_equal(attr(out, "perturbation_tests")$status, "not-computable")
    })
})

test_that("memoisation preserves each pathway outcome and outage status", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = simpleError("offline failure")), {
        x <- rbind(pathway_test_input(), pathway_test_input())
        x$Description <- c("first", "second")
        warnings <- character()
        out <- withCallingHandlers(pathway_evidence(x, drivers = "TF1",
            signed = pathway_test_signed), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        expect_true(any(grepl("All 1 GPSAdb calls failed", warnings)))
        expect_equal(attr(out, "network_calls"), 1L)
        expect_equal(attr(out, "perturbation_tests")$pathway, c("first", "second"))
        expect_equal(attr(out, "perturbation_tests")$status, rep("failed", 2))
        expect_equal(out$reason, rep("no-test-completed", 2))
    })
})

test_that("non-significance is not described as proof of falsity", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = data.frame(es = -1, padj = 0.5)), {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = pathway_test_signed))
        expect_equal(out$verdict, "coherent-refuted") # legacy API key
        expect_equal(out$reason, "not-significant")
        expect_equal(attr(out, "perturbation_tests")$status, "not-significant")
        expect_match(out$statement, "found no significant movement")
    })
})

test_that("requested candidates survive unresolved datasets and empty results", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    requested <- c("TF1", "UNRESOLVED")
    with_pathway_stubs(ds, list(D1 = data.frame(es = -1, padj = 0.001)), {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(),
            drivers = requested, signed = pathway_test_signed))
        expect_identical(attr(out, "drivers"), requested)
        empty <- suppressWarnings(pathway_evidence(pathway_test_input()[0, ], drivers = requested))
        expect_identical(attr(empty, "drivers"), requested)
        expect_null(attr(.pathway_evidence_empty(), "drivers"))
    })
})

test_that("row subsets scope the ledger without losing requested candidates", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = data.frame(es = -1, padj = 0.001)), {
        x <- rbind(pathway_test_input(), pathway_test_input())
        x$Description <- c("first", "second")
        out <- suppressWarnings(pathway_evidence(x, drivers = c("TF1", "UNRESOLVED"),
            signed = pathway_test_signed))
        expect_equal(attr(out[1, ], "perturbation_tests")$pathway, "first")
        expect_equal(attr(out[2, , drop = FALSE], "perturbation_tests")$pathway, "second")
        expect_equal(nrow(attr(out[FALSE, ], "perturbation_tests")), 0L)
        expect_identical(attr(out[FALSE, ], "drivers"), c("TF1", "UNRESOLVED"))
        expect_equal(attr(out[1, "level", drop = FALSE], "perturbation_tests")$pathway, "first")
        expect_equal(nrow(attr(out["pathway"], "perturbation_tests")), 2L)
    })
})

test_that("blank missing or padded cell lines cannot fabricate replication", {
    for (other in c("", " ", NA_character_, " cellA ")) {
        ds <- data.frame(driver = "TF1", dataset_id = c("D1", "D2"),
                         cell_line = c("cellA", other))
        got <- list(D1 = data.frame(es = -1, padj = 0.001),
                    D2 = data.frame(es = -1, padj = 0.002))
        with_pathway_stubs(ds, got, {
            warnings <- character()
            out <- withCallingHandlers(pathway_evidence(pathway_test_input(), drivers = "TF1",
                signed = pathway_test_signed), warning = function(w) {
                    warnings <<- c(warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                })
            expect_equal(out$level, "B")
            expect_equal(out$n_cell_lines_supported, 1L)
            expect_true(any(grepl("Level A is unreachable", warnings)))
        })
    }
})

test_that("conflicting identities for a dataset key are rejected", {
    ds <- data.frame(driver = "TF1", dataset_id = c("D1", "D1", "D2"),
                     cell_line = c("cellA", "cellB", "cellC"))
    got <- list(D1 = simpleError("ambiguous dataset must not run"),
                D2 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        warnings <- character()
        out <- withCallingHandlers(pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = pathway_test_signed), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        expect_true(any(grepl("conflicting cell-line identities: D1", warnings)))
        expect_equal(out$level, "B")
        expect_equal(out$n_cell_lines_supported, 1L)
        expect_equal(attr(out, "network_calls"), 1L)
        expect_equal(attr(out, "perturbation_tests")$dataset, "D2")
    })
})

test_that("Entrez conversion includes IDs beyond the fifth pathway", {
    skip_if_not_installed("org.Hs.eg.db")
    ns <- asNamespace("AnnotationDbi")
    original <- get("select", ns)
    on.exit(assignInNamespace("select", original, ns = "AnnotationDbi"), add = TRUE)
    seen <- NULL
    assignInNamespace("select", function(x, keys, ...) {
        seen <<- keys
        data.frame(ENTREZID = keys, SYMBOL = paste0("G", keys))
    }, ns = "AnnotationDbi")
    x <- data.frame(Description = paste0("path", 1:6), geneID = as.character(1:6))
    out <- .pathway_ids_to_symbols(x, c(G1 = 1))
    expect_identical(seen, as.character(1:6))
    expect_equal(nrow(out), 6L)
    expect_identical(out$geneID, paste0("G", 1:6))
})

test_that("dataset resolution accepts every removal method, not only 'knockdown'", {
    ## GPSAdb spells removal siRNA/shRNA/CRISPRi/knockout. A literal
    ## method = "knockdown" filter reported STAT1 -- CRISPRi plus knockout,
    ## no literal knockdown -- as having no data, the same false negative the
    ## inventory lookup was fixed for. Overexpression stays excluded.
    d <- data.frame(
        dataset_id = c("d1", "d2", "d3", "d4", "d5"),
        cell_line = c("A", "B", "C", "D", "E"),
        method = c("CRISPRi", "knockout", "siRNA", "overexpression", "shRNA"),
        in_reference = c(TRUE, TRUE, TRUE, TRUE, FALSE),
        n_samples = c(3, 2, 4, 9, 1), stringsAsFactors = FALSE)
    got <- .pathway_usable_datasets(d)
    expect_identical(got$dataset_id, c("d1", "d2", "d3"))
    expect_false(any(grepl("overexpress", got$method)))
})

test_that("unusable dataset tables are empty, not fabricated data", {
    expect_null(.pathway_usable_datasets(NULL))
    expect_null(.pathway_usable_datasets(data.frame()))
    expect_null(.pathway_usable_datasets(data.frame(dataset_id = "d1")))
    d <- data.frame(dataset_id = c("d1", "d1"), cell_line = c("A", " "),
                    method = "siRNA", in_reference = TRUE, stringsAsFactors = FALSE)
    expect_null(suppressWarnings(.pathway_usable_datasets(d)))
    blank <- data.frame(dataset_id = c("d1", "d2"), cell_line = c("A", "  "),
                        method = "siRNA", in_reference = TRUE, stringsAsFactors = FALSE)
    out <- .pathway_usable_datasets(blank)
    expect_equal(nrow(out), 2L)
    expect_true(is.na(out$cell_line[2]))
})

test_that("a failed dataset lookup is not reported as absent data", {
    ns <- asNamespace("gpsadb")
    original <- get("gpsadb_search_datasets", ns)
    on.exit(assignInNamespace("gpsadb_search_datasets", original, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_search_datasets",
                      function(...) stop("simulated outage"), ns = "gpsadb")
    warnings <- character()
    out <- withCallingHandlers(
        pathway_evidence(pathway_test_input(), drivers = "TF1",
                         signed = pathway_test_signed),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_equal(out$reason, "lookup-failed")
    expect_match(out$statement, "lookup failure, not evidence")
    expect_true(any(grepl("lookup failure", warnings)))
    expect_false(any(grepl("No usable perturbation data found", warnings)))
    lookup <- attr(out, "dataset_lookup")
    expect_equal(lookup$driver, "TF1")
    expect_equal(lookup$status, "failed")
})

test_that("an empty lookup stays distinguishable from a failed one", {
    ns <- asNamespace("gpsadb")
    original <- get("gpsadb_search_datasets", ns)
    on.exit(assignInNamespace("gpsadb_search_datasets", original, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_search_datasets",
                      function(...) list(data.frame(dataset_id = character(),
                          cell_line = character(), method = character(),
                          in_reference = logical(), n_samples = integer())),
                      ns = "gpsadb")
    warnings <- character()
    out <- withCallingHandlers(
        pathway_evidence(pathway_test_input(), drivers = "TF1",
                         signed = pathway_test_signed),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_equal(out$reason, "no-dataset")
    expect_true(any(grepl("No usable perturbation data found", warnings)))
    expect_false(any(grepl("lookup failure", warnings)))
    expect_equal(attr(out, "dataset_lookup")$status, "empty")
})

test_that("a partial lookup failure is disclosed even when another driver resolves", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    ns <- asNamespace("gpsadb")
    original <- get("gpsadb_search_datasets", ns)
    on.exit(assignInNamespace("gpsadb_search_datasets", original, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        if (identical(perturbed_gene, "TF2")) stop("simulated outage for TF2")
        z <- ds[ds$driver == perturbed_gene, , drop = FALSE]
        z$n_samples <- 4L
        z$in_reference <- TRUE
        list(z)
    }, ns = "gpsadb")
    original_test <- get("gpsadb_test_gene_set", ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", original_test, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) list(got[[id]]),
                      ns = "gpsadb")
    warnings <- character()
    out <- withCallingHandlers(
        pathway_evidence(pathway_test_input(), drivers = c("TF1", "TF2"),
                         signed = pathway_test_signed),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_true(any(grepl("could not be queried for `TF2`", warnings)))
    expect_true(any(grepl("lookup failure", warnings)))
    ## The single-cell-line warning must name the resolved candidate only: TF2
    ## was never queried, so it cannot have evidence in one cell line.
    single <- warnings[grepl("single cell line", warnings)]
    expect_true(length(single) >= 1L)
    expect_true(any(grepl("TF1", single)))
    expect_false(any(grepl("TF2", single)))
    lookup <- attr(out, "dataset_lookup")
    expect_equal(lookup$status[lookup$driver == "TF1"], "available")
    expect_equal(lookup$status[lookup$driver == "TF2"], "failed")
})

test_that("missing fold changes are not blamed on the identifier space", {
    ## A documented call: direction from pathway_direction/NES, no `signed`.
    ## The perturbation tests run and support the direction, but member support
    ## is unknown, so the row cannot pass the support floor. It must say the
    ## fold changes are missing, not that the members are.
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1"))
        expect_equal(out$level, "C")
        expect_equal(out$reason, "no-fold-changes")
        expect_match(out$statement, "No signed fold changes")
        expect_false(grepl("signed gene list", out$statement))
        expect_true(is.na(out$support))
        expect_equal(attr(out, "perturbation_tests")$status, "supports")
    })
})

test_that("a column-only subset keeps the retained ledger", {
    two <- data.frame(Description = c("p1", "p2"), geneID = "A/B/C/D/E",
                      pathway_direction = "up", stringsAsFactors = FALSE)
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        pe <- suppressWarnings(pathway_evidence(two, drivers = "TF1",
            signed = pathway_test_signed))
        reordered <- pe[, rev(names(pe))]
        expect_s3_class(reordered, "pathway_evidence")
        expect_equal(nrow(attr(reordered, "perturbation_tests")), 2L)
        expect_identical(attr(reordered, "perturbation_tests"),
                         attr(pe, "perturbation_tests"))
        rep <- mechanism_report(reordered, inventory = data.frame(perturbed_gene = "TF1"))
        expect_equal(nrow(rep$perturbation_tests), 2L)
        expect_equal(nrow(rep$candidate_evidence), 2L)
        ## Dropping a required column still falls back to a plain data.frame.
        slim <- pe[, c("pathway", "level")]
        expect_false(inherits(slim, "pathway_evidence"))
    })
})
