## 发现门（mechanism_discover）的判定与守卫。
##
## 这一层不联网、不依赖外部表，所以全部可测。

.make_net <- function(df) data.frame(source = df[[1]], target = df[[2]],
                                     weight = df[[3]], stringsAsFactors = FALSE)

test_that("a factor whose targets moved ranks first", {
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(c(rep(1.2, 6), rep(0, 34)), genes)
    ## 显式写 data.frame：`.make_net(list(...))` 会因为向量长度不等而被回收，
    ## 把两个 TF 的靶基因交错成同一批（测试第一版就是这么错的）。
    net <- data.frame(source = rep(c("TF1", "TF2"), each = 6),
                      target = c(genes[1:6], genes[7:12]), weight = 1,
                      stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 2)
    expect_equal(out$source[1], "TF1")
    expect_equal(out$direction[1], "activates")
    expect_gt(out$score[1], out$score[2])
})

test_that("self-targets are dropped, and that is what keeps knockouts honest", {
    ## 敲除实验里 TF 自己的 log2FC 最强；自环会把它的分数顶到第一。
    ## 这不是理论问题：CollecTRI 里 62% 的 TF 带自环，删掉值 5 个点。
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(c(-3, rep(1.2, 6), rep(0, 33)), genes)
    ## 自环在真实数据里是"TF 把自己列为靶"——所以 TF 名必须就是基因名之一。
    ## g1 是 TF1，它的靶 = 自己(g1) + g2..g7；TF2 的靶 = g2..g7。
    net <- data.frame(
        source = c("g1", rep("g1", 6), rep("g2", 6)),
        target = c("g1", genes[2:7], genes[2:7]),
        weight = 1, stringsAsFactors = FALSE)
    expect_message(out <- mechanism_discover(signed, net, top = 2), "self-target")
    ## 去掉自环后 TF1 与 TF2 靶集相同、分数并列；关键是 TF1 不再靠自身 log2FC 领先
    expect_equal(out$n_targets[out$source == "g1"], 6L)
    keep <- suppressMessages(mechanism_discover(signed, net, top = 2, keep_self_targets = TRUE))
    expect_equal(keep$n_targets[keep$source == "g1"], 7L)
})

test_that("weights carry the direction of regulation", {
    genes <- paste0("g", 1:40)
    ## 靶基因整体下调，但权重是"抑制靶"（-1）→ 该因子应当判为激活
    signed <- stats::setNames(c(rep(-1.2, 6), rep(0, 34)), genes)
    net <- .make_net(list("TF1", genes[1:6], rep(-1, 6)))
    out <- mechanism_discover(signed, net, top = 1)
    expect_equal(out$direction[1], "activates")
    expect_gt(out$score[1], 0)
})

test_that("factors with too few measured targets are skipped", {
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(stats::rnorm(40), genes)
    net <- data.frame(source = c(rep("TF1", 3), rep("TF2", 17)),
                      target = c(genes[1:3], genes[4:20]), weight = 1,
                      stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 5, min_targets = 5)
    expect_equal(nrow(out), 1)
    expect_equal(out$source, "TF2")
})

test_that("targets outside the measured genes are dropped, not zero-filled", {
    ## 零填会把"没测到"当成"没变化"，分数被稀释。必须按测到的交集算。
    genes <- paste0("g", 1:40)
    ## 只有 5 个基因被网络覆盖，所以它们之间必须有差异，否则方差为 0（这是正确的守卫）
    signed <- stats::setNames(c(1, 1.1, 0.9, 1.05, 0.95, rep(0, 35)), genes)
    net <- data.frame(source = "TF1",
                      target = c(genes[1:5], "not_measured_1", "not_measured_2"),
                      weight = 1, stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 1)
    expect_equal(out$n_targets, 5)
    expect_equal(attr(out, "n_genes"), 40)
})

test_that("a network with no overlap fails loudly", {
    signed <- stats::setNames(c(a = 1, b = 2), c("a", "b"))
    net <- .make_net(list("TF1", c("x", "y", "z"), rep(1, 3)))
    expect_error(mechanism_discover(signed, net), "None of the network targets")
})

test_that("the shortlist is capped and ordered", {
    genes <- paste0("g", 1:60)
    signed <- stats::setNames(stats::rnorm(60), genes)
    net <- .make_net(list(rep(paste0("TF", 1:5), each = 8),
                          rep(genes[1:40], 5), rep(1, 40)))
    out <- mechanism_discover(signed, net, top = 3)
    expect_equal(nrow(out), 3)
    expect_equal(out$rank, 1:3)
    expect_true(all(diff(out$score) <= 0))
})
