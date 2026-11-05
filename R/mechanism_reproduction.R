## 复现率：一条单细胞系的 B 级结论，到底有没有第二份数据撑住？
##
## 为什么要这一层。11 个已知真值的数据集里，停在 B 级（`single-cell-line`）的两行
## ——airway 与 GSE286890 的 GR——**都指向正确的驱动者**。所以 B 不等于错，B 是
## "证据只有一份"。但"只有一份"与"多份都指向同一方向"在读法上完全不同，而原来的
## 输出对两者用同一句措辞。
##
## 判据来自 `mechanism-benchmark-spec.md` §4：**一致性**（concordance），不是
## "与通路方向一致"。后者是符号错误——激活型驱动者被敲低时，通路应当**下调**，
## 与通路方向**相反**。那会把真驱动者算成 0 复现，而且越真的越接近 0
## （实测 HIF1A × 缺氧通路：旧口径 0.00，新口径 1.00，majority_sign = −1）。
##
## 方向判定不经过 `pathway_evidence()` 自己的 es：直接用扰动谱的成员均值符号，
## 否则是拿同一统计量校准自己的措辞（第一期 §7.10 记过一次自引用）。

#' Concordance of a pathway shift across independent perturbation datasets
#'
#' For each (driver, pathway) pair, takes the mean log2 fold change of the
#' pathway members in every available perturbation dataset of that driver, and
#' reports how many of them move the pathway the same way.
#'
#' This is the input [pathway_evidence()] needs to decide how strongly a
#' single-cell-line level-B row may be worded.
#'
#' @param profiles Numeric matrix, rows are genes (symbols), columns are
#'   perturbation dataset identifiers.
#' @param meta Data frame with columns `dataset` and `driver`. Optional
#'   `biosample` (for the independent-sample count) and `logFC` (the knockdown
#'   strength, used by `max_log2fc`).
#' @param pathways Named list: pathway name to a character vector of member
#'   genes. [pathway_evidence()] builds this from the enrichment result.
#' @param drivers Character vector of drivers to score. Defaults to every
#'   driver present in `meta`.
#' @param min_datasets Integer. Fewer datasets than this cannot show concordance
#'   (one dataset is trivially concordant), so the pair is reported as
#'   `not-computable`.
#' @param min_biosamples Integer. A row only counts as reproduced when the
#'   supporting datasets come from at least this many distinct biosamples.
#' @param min_genes Integer. A pathway with fewer measured members than this is
#'   reported as `not-computable`: the mean of one or two genes is noise, not a
#'   shift. The etoposide dataset in our own benchmark has pathways of 1-2
#'   members, which is exactly why this floor exists.
#' @param max_log2fc Optional numeric. When given and `meta` has a `logFC`
#'   column, datasets weaker than this knockdown are dropped before scoring.
#'
#' @return A data frame with one row per (driver, pathway): `driver`, `pathway`,
#'   `k` (datasets used), `n_biosamples`, `n_genes` (members found in
#'   `profiles`), `concordance`, `majority_sign`, and `grade` which is one of
#'   `reproduces`, `consistent`, `not-reproduced` or `not-computable`.
#'
#' @examples
#' set.seed(1)
#' profiles <- matrix(rnorm(60), nrow = 20,
#'                    dimnames = list(paste0("g", 1:20), c("d1", "d2", "d3")))
#' profiles[1:5, ] <- -1
#' meta <- data.frame(dataset = c("d1", "d2", "d3"),
#'                    driver = "TF1", biosample = c("A", "B", "B"))
#' mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:5)))
#'
#' @export
mechanism_reproduction <- function(profiles, meta, pathways, drivers = NULL,
                                   min_datasets = 2L, min_biosamples = 2L,
                                   min_genes = 3L, max_log2fc = NULL) {
    if (!is.matrix(profiles) || !is.numeric(profiles))
        stop("`profiles` must be a numeric matrix.", call. = FALSE)
    if (is.null(rownames(profiles)) || is.null(colnames(profiles)))
        stop("`profiles` must have gene row names and dataset column names.",
             call. = FALSE)
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (!all(c("dataset", "driver") %in% names(meta)))
        stop("`meta` must have `dataset` and `driver` columns.", call. = FALSE)
    if (!is.list(pathways) || is.null(names(pathways)))
        stop("`pathways` must be a named list of gene vectors.", call. = FALSE)
    meta$dataset <- as.character(meta$dataset)
    meta$driver <- as.character(meta$driver)
    meta <- meta[meta$dataset %in% colnames(profiles), , drop = FALSE]
    ## 敲低强度过滤。放在取子集之前，否则 k 会把弱敲低的实验算进去——
    ## 实测 12 个驱动者里有 3 个（ARNT / ESR1 / NR3C1）全部实验都弱于 −1，
    ## 不过滤时它们各得 4 / 6 / 1 份，过滤后归零。这个差别必须由调用方看见，
    ## 所以 max_log2fc 是参数而不是硬编码。
    if (!is.null(max_log2fc) && "logFC" %in% names(meta)) {
        keep <- !is.na(meta$logFC) & meta$logFC <= max_log2fc
        meta <- meta[keep, , drop = FALSE]
    }
    if (is.null(drivers)) drivers <- unique(meta$driver)
    drivers <- as.character(drivers)
    has_bio <- "biosample" %in% names(meta)

    rows <- vector("list", length(pathways) * length(drivers))
    n <- 0L
    for (p in names(pathways)) {
        members <- intersect(as.character(pathways[[p]]), rownames(profiles))
        for (dr in drivers) {
            n <- n + 1L
            sel <- meta[meta$driver == dr, , drop = FALSE]
            k <- nrow(sel)
            nb <- if (has_bio) length(unique(sel$biosample)) else NA_integer_
            if (k < min_datasets || length(members) < min_genes) {
                rows[[n]] <- data.frame(
                    driver = dr, pathway = p, k = k, n_biosamples = nb,
                    n_genes = length(members), concordance = NA_real_,
                    majority_sign = NA_integer_, grade = "not-computable",
                    stringsAsFactors = FALSE)
                next
            }
            sub <- profiles[members, sel$dataset, drop = FALSE]
            shift <- colMeans(sub, na.rm = TRUE)
            shift <- shift[is.finite(shift)]
            maj <- if (!length(shift)) NA_integer_ else sign(sum(shift))
            if (is.na(maj) || maj == 0) {
                conc <- NA_real_
                grade <- "not-computable"
                maj <- NA_integer_
            } else {
                conc <- mean(sign(shift) == maj)
                grade <- if (conc >= 0.9 && k >= 3L &&
                             (is.na(nb) || nb >= min_biosamples)) "reproduces"
                         else if (conc >= 0.5) "consistent"
                         else "not-reproduced"
            }
            rows[[n]] <- data.frame(
                driver = dr, pathway = p, k = k, n_biosamples = nb,
                n_genes = length(members), concordance = conc,
                majority_sign = maj, grade = grade, stringsAsFactors = FALSE)
        }
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    attr(out, "min_datasets") <- as.integer(min_datasets)
    out
}
