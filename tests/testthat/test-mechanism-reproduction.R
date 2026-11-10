## 复现率（mechanism_reproduction）与它在 pathway_evidence 里的用法。
##
## 这一层的测试只关心**判定规则**，不碰网络：输入是构造好的矩阵，
## 输出是可预期的分级。真实数据的复现率由 `scripts/` 里的离线脚本算。

test_that("concordance is measured against the majority, not against the pathway", {
    ## 激活型驱动者被敲低 → 成员应当**下调**。若拿"与通路方向一致"当判据，
    ## 真驱动者会得到 0 复现（实跑 HIF1A：旧口径 0.00、新口径 1.00）。
    profiles <- matrix(0, nrow = 6, ncol = 3,
                       dimnames = list(paste0("g", 1:6), c("d1", "d2", "d3")))
    profiles[1:4, ] <- -1        # 四分之三的成员下调，三份数据一致
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)))
    expect_equal(out$concordance, 1)
    expect_equal(out$majority_sign, -1L)
    expect_equal(out$grade, "reproduces")
})

test_that("a single dataset is not-computable rather than perfectly concordant", {
    ## k = 1 时一致性恒为 1。若不特判，"只有一份数据"会被算成"完美复现"，
    ## 恰好把最该保守的情形说成最强的结论。
    profiles <- matrix(-1, nrow = 4, ncol = 1,
                       dimnames = list(paste0("g", 1:4), "d1"))
    meta <- data.frame(dataset = "d1", driver = "TF1", biosample = "A")
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)))
    expect_true(is.na(out$concordance))
    expect_equal(out$grade, "not-computable")
})

test_that("runs that disagree are graded not-reproduced", {
    profiles <- matrix(c(-1, -1, 1), nrow = 1, ncol = 3,
                       dimnames = list("g1", c("d1", "d2", "d3")))
    profiles <- rbind(profiles, g2 = c(-1, 1, 1), g3 = c(-1, 1, 1))
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(profiles, meta, list(pw = c("g1", "g2", "g3")))
    ## 成员均值：d1 = -1, d2 = +1/3, d3 = +1/3 → 多数派为正，一致性 2/3
    expect_equal(out$concordance, 2 / 3)
    expect_equal(out$grade, "consistent")
})

test_that("max_log2fc drops weak knockdowns before counting datasets", {
    ## 实测：ARNT / ESR1 / NR3C1 在 KnockTF 里的全部实验都弱于 -1，
    ## 不过滤时各得 4 / 6 / 1 份，过滤后归零。这个差别必须能被调用方看见。
    profiles <- matrix(-1, nrow = 4, ncol = 3,
                       dimnames = list(paste0("g", 1:4), c("d1", "d2", "d3")))
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"),
                       logFC = c(-2, -0.2, -0.3))
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)),
                                  max_log2fc = -1)
    expect_equal(out$k, 1L)
    expect_equal(out$grade, "not-computable")
})

test_that("reproduces needs three datasets from at least two biosamples", {
    profiles <- matrix(-1, nrow = 4, ncol = 3,
                       dimnames = list(paste0("g", 1:4), c("d1", "d2", "d3")))
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "A", "A"))
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)))
    expect_equal(out$concordance, 1)
    ## 三份同细胞系不算独立复现
    expect_equal(out$grade, "consistent")
})

test_that("the level-B wording has three distinguishable branches", {
    ## 不联网：直接测措辞辅助函数。三个分支必须彼此可区分——
    ## 若"无第二份数据"与"不复现"共用一句，读者就分不清"没人验过"与"验过但失败了"。
    base <- data.frame(driver = "TF1", pathway = "pw", k = 2L, n_biosamples = 2L,
                       n_genes = 3L, concordance = 1, majority_sign = -1L,
                       grade = "reproduces", stringsAsFactors = FALSE)
    unrep <- base; unrep$concordance <- 0; unrep$grade <- "not-reproduced"
    nocov <- base; nocov$k <- 1L; nocov$concordance <- NA_real_
    nocov$grade <- "not-computable"

    st_rep <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, base)
    expect_match(st_rep, "reproduces across 2 independent datasets", fixed = TRUE)
    expect_match(st_rep, "2/2 agree", fixed = TRUE)

    st_unrep <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, unrep)
    expect_match(st_unrep, "unverified", fixed = TRUE)
    expect_match(st_unrep, "0/2 agree", fixed = TRUE)

    st_nocov <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, nocov)
    expect_match(st_nocov, "No independent reproduction is available", fixed = TRUE)

    ## 没有 reproduction 表时退回原措辞
    st_none <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, NULL)
    expect_match(st_none, "limit of the available data", fixed = TRUE)

    expect_false(identical(st_rep, st_unrep))
    expect_false(identical(st_unrep, st_nocov))
})

test_that("the majority sign is the mode of signs, not the sign of the sum", {
    ## 一份效应量特别大的数据不能翻掉两份较小但同向的。
    ## 旧实现用 sign(sum(shift))：sum(c(-0.1,-0.1,+5)) > 0，于是把
    ## "2/3 同向"判成 not-reproduced。
    prof <- matrix(0, nrow = 3, ncol = 3,
                   dimnames = list(c("g1", "g2", "g3"), c("d1", "d2", "d3")))
    prof[1, ] <- c(-0.1, -0.1, 5)
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(prof, meta, list(pw = c("g1", "g2", "g3")))
    expect_equal(out$majority_sign, -1L)
    expect_equal(out$concordance, 2 / 3)
    expect_equal(out$grade, "consistent")
})

test_that("datasets with no measured members do not inflate k", {
    prof <- matrix(0, nrow = 3, ncol = 3,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c")))
    prof[1, ] <- c(-1, -1, -1)
    prof[, 3] <- NA_real_
    meta <- data.frame(dataset = c("a", "b", "c"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(prof, meta, list(pw = c("g1", "g2", "g3")))
    expect_equal(out$k, 2L)
    expect_equal(out$n_biosamples, 2L)
    expect_equal(out$grade, "consistent")   # k < 3，不能判 reproduces
})

test_that("finite member floors apply separately to every dataset", {
    prof <- matrix(c(-1, NA, Inf), 3, 3,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c")))
    meta <- data.frame(dataset = c("a", "b", "c"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
    expect_equal(out$k, 0L)
    expect_equal(out$n_biosamples, 0L)
    expect_equal(out$grade, "not-computable")
    prof[, 1] <- -1
    out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
    expect_equal(out$k, 1L)
    expect_equal(out$grade, "not-computable")
})

test_that("missing biosamples never demonstrate independent reproduction", {
    prof <- matrix(-1, 3, 3,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c")))
    meta <- data.frame(dataset = c("a", "b", "c"), driver = "TF1")
    out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
    expect_true(is.na(out$n_biosamples))
    expect_equal(out$grade, "consistent")
    for (bio in list(c("A", NA, NA), c("A", "", "  "), rep(NA_character_, 3))) {
        meta$biosample <- bio
        out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
        expect_lte(out$n_biosamples, 1L)
        expect_equal(out$grade, "consistent")
    }
})

test_that("dataset identities cannot inflate reproduction counts", {
    prof <- matrix(-1, 3, 2,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b")))
    meta <- data.frame(dataset = c("a", "b"), driver = "TF1", biosample = c("A", "B"))
    out <- mechanism_reproduction(prof, rbind(meta, meta[1, ]), list(pw = rownames(prof)))
    expect_equal(out$k, 2L)
    expect_equal(out$grade, "consistent")
    conflict <- meta[1, ]; conflict$biosample <- "C"
    expect_error(mechanism_reproduction(prof, rbind(meta, conflict), list(pw = rownames(prof))),
                 "inconsistent duplicate")
    colnames(prof) <- c("a", "a")
    expect_error(mechanism_reproduction(prof, meta, list(pw = rownames(prof))), "unique")
})

test_that("reproduction rows are matched through driver aliases", {
    ## GPSAdb 返回 NR3C1，用户的表可能用别名 GR 建成。不解析会静默退回默认措辞，
    ## 而这恰好是 P4 里两个 B 级案例（airway / GSE286890）的驱动者。
    rep <- data.frame(driver = "GR", pathway = "pw", k = 3L,
                      n_biosamples = 3L, n_genes = 5L, concordance = 1,
                      majority_sign = -1L, grade = "reproduces",
                      stringsAsFactors = FALSE)
    hit <- .pe_reproduction_row("NR3C1", "pw", rep)
    expect_equal(nrow(hit), 1L)
    st <- .pe_single_line_statement("NR3C1", 5L, "D1", 0.01, hit)
    expect_match(st, "reproduces across 3 independent datasets", fixed = TRUE)
    ## 反方向也要通
    rep2 <- rep; rep2$driver <- "NR3C1"
    expect_equal(nrow(.pe_reproduction_row("GR", "pw", rep2)), 1L)
    expect_error(.pe_reproduction_row("NR3C1", "pw", data.frame(x = 1)),
                 "driver")
})
