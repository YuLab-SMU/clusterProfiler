## 通路模式的第一批测试。
##
## 这些测的是**判断逻辑**，不是网络。`pathway_evidence()` 在有驱动者时会调 GPSAdb，
## 所以涉及网络的用例一律 skip；离线用例覆盖方向推导、成员一致性、分级与陈述。

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
    expect_true("n_datasets_tested" %in% names(out))
    ## 若有显著支持但只有单一细胞系，不能是 A。
    if (out$n_datasets_tested > 0 && out$padj < 0.05 && out$n_cell_lines_supported < 2) {
        expect_false(out$level == "A")
        expect_equal(out$verdict, "single-dataset")
    }
})
