## 预检门的测试。
##
## 重点不是"能不能打印"，而是**每一个挡人的理由都对应一个可复现的输入**，
## 以及**不该挡的时候不挡**——一个总说"不行"的门和没有门一样没用。

.make_preflight_ref <- function(n_cell = 2L, n_prog = 20L) {
    ## 每个基因两份数据集（两份不同细胞系），这样"驱动者有 ≥2 个细胞系"这条
    ## 在健康输入上成立；n_cell = 1 时把两份都放到同一个细胞系上，用来触发那一条。
    g <- sprintf("TF%02d", seq_len(n_prog))
    meta <- data.frame(
        dataset_id = paste0("D", rep(seq_len(n_prog), each = 2), c("a", "b")),
        perturbed_gene = rep(g, each = 2),
        cell_line = if (n_cell >= 2L) rep(c("A", "B"), n_prog) else rep("A", 2 * n_prog),
        tissue = "T", method = "siRNA", n_samples = 4L,
        stringsAsFactors = FALSE)
    genes <- sprintf("G%03d", seq_len(100))
    structure(list(scores = matrix(0, length(genes), nrow(meta),
                                   dimnames = list(genes, meta$dataset_id)),
                   meta = meta, coverage = NULL, overall_coverage = 1,
                   global = stats::setNames(rep(0, length(genes)), genes),
                   global_var_pct = 1, global_imputed = FALSE,
                   loadings = data.frame(gene = genes[1], loading = 0),
                   genes = genes, failed = character(0)),
              class = "mechanismReference")
}

.make_ego <- function(n = 5L, members = 20L) {
    data.frame(Description = paste("p", seq_len(n)),
               ID = paste0("GO:", seq_len(n)),
               geneID = rep(paste(sprintf("G%03d", seq_len(members)), collapse = "/"), n),
               stringsAsFactors = FALSE)
}

test_that("no pathways means no tier above none", {
    empty <- data.frame(Description = character(), ID = character(),
                        geneID = character(), stringsAsFactors = FALSE)
    pf <- mechanism_preflight(empty)
    expect_equal(pf$tier, "none")
    expect_match(pf$reasons[1], "no-pathways")
    expect_output(print(pf), "none")
})

test_that("too few up genes stops attribution", {
    ego <- .make_ego()
    signed <- stats::setNames(c(rep(2, 20), rep(0.1, 80)), sprintf("G%03d", seq_len(100)))
    pf <- mechanism_preflight(ego, signed = signed)
    expect_equal(pf$tier, "description")
    expect_match(paste(pf$reasons, collapse = " "), "few-up-genes")
    expect_equal(pf$metrics$n_up, 20L)
})

test_that("a driver with one cell line stops attribution", {
    ego <- .make_ego()
    signed <- stats::setNames(rep(2, 100), sprintf("G%03d", seq_len(100)))
    ref <- .make_preflight_ref(n_cell = 1L)
    pf <- mechanism_preflight(ego, signed = signed, reference = ref, drivers = "TF01")
    expect_equal(pf$tier, "description")
    expect_match(paste(pf$reasons, collapse = " "), "single-cell-line")
})

test_that("small pathways stop attribution even with many genes", {
    ego <- .make_ego(members = 3L)
    signed <- stats::setNames(rep(2, 100), sprintf("G%03d", seq_len(100)))
    pf <- mechanism_preflight(ego, signed = signed)
    expect_equal(pf$tier, "description")
    expect_match(paste(pf$reasons, collapse = " "), "small-pathways")
})

test_that("a thin reference stops attribution", {
    ego <- .make_ego()
    signed <- stats::setNames(rep(2, 100), sprintf("G%03d", seq_len(100)))
    ref <- .make_preflight_ref(n_prog = 4L)
    pf <- mechanism_preflight(ego, signed = signed, reference = ref)
    expect_equal(pf$tier, "description")
    expect_match(paste(pf$reasons, collapse = " "), "thin-reference")
})

test_that("a healthy input is not blocked", {
    ego <- .make_ego()
    signed <- stats::setNames(rep(2, 100), sprintf("G%03d", seq_len(100)))
    ref <- .make_preflight_ref()
    pf <- mechanism_preflight(ego, signed = signed, reference = ref, drivers = "TF01")
    expect_equal(pf$tier, "attribution")
    expect_length(pf$reasons, 0L)
})

test_that("the tier is the worst reason, not the first", {
    ego <- .make_ego(members = 3L)
    signed <- stats::setNames(c(rep(2, 10), rep(0.1, 90)), sprintf("G%03d", seq_len(100)))
    pf <- mechanism_preflight(ego, signed = signed)
    expect_equal(pf$tier, "description")
    expect_true(length(pf$reasons) >= 2L)
    expect_match(paste(pf$reasons, collapse = " "), "few-up-genes")
    expect_match(paste(pf$reasons, collapse = " "), "small-pathways")
})

test_that("one under-covered candidate does not block the whole run", {
    ## 4 个候选里 3 个有 2 个细胞系、1 个只有 1 个：整批仍然可以归属，
    ## 只是那一个候选自己不该被认领。取 min 会把整批都挡掉。
    ref <- .make_preflight_ref()
    ref$meta$cell_line[ref$meta$perturbed_gene == "TF20"] <- "A"
    ego <- .make_ego()
    signed <- stats::setNames(rep(2, 100), sprintf("G%03d", seq_len(100)))
    pf <- mechanism_preflight(ego, signed = signed, reference = ref,
                              drivers = c("TF01", "TF02", "TF20"))
    expect_equal(pf$tier, "attribution")
    expect_equal(pf$metrics$min_cell_lines, 1L)
    expect_equal(pf$metrics$n_cell_lines, 2L)
})
