## 预检门：在跑任何归属之前，先判断这份数据支不支持归属。
##
## 为什么要有这一层。十一个数据集里五个失败，而失败原因**在跑之前就能看出来**：
## 两个数据集 GO 一条都没有、一个只有 23 个上调基因、一个驱动者在库里只有一份
## 数据集。原来的做法是照跑一遍，然后输出一排 `C unsupported` ——读者分不清
## "没找到"与"没法找"。
##
## 判据来自 `pathway-mechanism-design.md` §7.9 已经记录的适用边界（差异基因要够多、
## 通路成员要够、驱动者要有 ≥2 个细胞系），**不是**从这次对比的结果反推的。
## 但那也意味着：§7.9 只给了定性说法（"23 明显不够、351 明显够"），下面这三个
## 阈值是**在两个已知点之间取的数，没有测过**。所以门输出的是**档位 + 依据的数字**，
## 而不是一个"合格"的判决；调用方可以按自己的数据调阈值。

## 上调基因数：用 log2FC > 1 计数。
##
## 严格的写法要 padj < 0.05 且 |log2FC| >= 1，但那需要调用方传 padj，而这一层
## 拿到的只有 `signed`。用 log2FC > 1 会把边界放宽一点（缺氧：351 严格 → 约 400），
## 对"够不够"这个判断没有影响，但要在文档里说清不是同一个数。
.mechanism_preflight_counts <- function(x, signed) {
    d <- if (is.data.frame(x)) x else as.data.frame(x)
    ## 优先用富集分析**真正吃进去**的那批基因（`enrichResult` 的 `@gene`）。
    ## 退回 `signed > 1` 会把数目放大好几倍（缺氧：严格 351，退回口径 881），
    ## 于是"上调基因太少"这条判据几乎永不触发。
    n_up <- if (inherits(x, "enrichResult") && length(methods::slot(x, "gene"))) {
        length(methods::slot(x, "gene"))
    } else if (!is.null(signed) && length(signed)) {
        sum(signed > 1, na.rm = TRUE)
    } else NA_integer_
    list(n_pathways = nrow(d), n_up = n_up, d = d)
}

#' Decide whether a dataset can support driver attribution at all
#'
#' Run this before [mechanism_decompose()] or [pathway_evidence()]. It returns a
#' tier rather than a verdict: `attribution` when the data can carry a candidate
#' driver claim, `description` when it can only describe what moved, and `none`
#' when even that is not possible.
#'
#' The boundaries are the ones already written down in the pathway-mechanism
#' design (differential genes, pathway size, driver coverage), not thresholds
#' fitted to a benchmark. They are still crude: the design notes record that 23
#' up-regulated genes are clearly too few and 351 clearly enough, so the cut
#' points below sit between two measured points rather than on a measured one.
#' The numbers behind the tier are returned with it, so a caller can move the
#' cut points and see which datasets change tier.
#'
#' @param x An enrichment result, or a data.frame of pathways.
#' @param signed Named numeric of log fold changes for every gene measured.
#' @param reference Optional `mechanismReference`; needed to judge whether the
#'   pathways have enough covered members and whether the candidates have data in
#'   more than one cell line.
#' @param drivers Optional candidate drivers. When given, a candidate with data
#'   in a single cell line caps the whole run at `description`, because a driver
#'   that cannot be seen in a second context cannot be separated from an artefact
#'   of the first.
#' @param min_up Up-regulated genes required before attribution is offered.
#' @param min_members Median pathway size required. This is the pathway's own
#'   size, not how many of its members the reference happens to cover: the
#'   12-member floor recorded in the design belongs to the perturbation test's
#'   p-value, and applying it here blocked three datasets that the decomposition
#'   handled fine.
#' @param min_cell_lines Cell lines a driver needs before attribution is offered.
#' @param min_background Comparable reference programmes required for the
#'   specificity percentile to mean anything.
#' @return A `mechanismPreflight` list with `tier`, `reasons`, `metrics` and
#'   `counts`.
#' @seealso [mechanism_decompose()], [mechanism_report()]
#' @export
mechanism_preflight <- function(x, signed = NULL, reference = NULL, drivers = NULL,
                                min_up = 50L, min_members = 5L,
                                min_cell_lines = 2L, min_background = 10L) {
    cnt <- .mechanism_preflight_counts(x, signed)
    reasons <- character(0)
    if (!cnt$n_pathways) {
        return(.mechanism_preflight_out(
            "none", c("no-pathways: the enrichment result holds no pathways, so there is nothing to describe or attribute."),
            cnt, NA_real_, NA_integer_, NA_integer_))
    }

    ## 通路大小按**通路自身**算，不按"覆盖成员数"。
    ##
    ## 第一版用的是覆盖成员数、阈值 8，结果把三个成功的数据集也挡了（它们的通路
    ## 中位覆盖只有 5–7）。那个 8 是从 §7.9 抄来的，而 §7.9 的"成员数要 12 以上"
    ## 说的是 `test_gene_set` 算不出 p 值——**那是旧统计量的边界，不是这条的**。
    ## 形状相关在 5 个成员的通路上照样能算。把旧边界搬到新统计量上，是类别错误。
    sz <- vapply(seq_len(nrow(cnt$d)), function(i) length(.pathway_members(cnt$d[i, ])), 0L)
    med_members <- stats::median(sz, na.rm = TRUE)

    n_cell_lines <- NA_integer_
    if (length(drivers)) {
        ## 有参考矩阵就从矩阵里数；没有就现场查一次库存（每个驱动者一次调用）。
        ## 没有这一步，旧实现的用户（不建矩阵）永远拿不到"单细胞系"这条判据——
        ## 而它恰恰是 11 个数据集里两个失败的唯一原因。
        per <- if (!is.null(reference) && inherits(reference, "mechanismReference")) {
            meta <- reference$meta
            prog <- stats::setNames(as.character(meta$perturbed_gene), meta$dataset_id)
            ## 别名要落到真有数据的基因上：`NFE2L2` 的别名 `NRF2` 没有数据集，
            ## 不剔除会得到 0 个细胞系，把舒拉芬误判成单细胞系。
            present <- intersect(drivers, unique(prog))
            vapply(present, function(g) {
                ids <- names(prog)[prog %in% g]
                length(unique(stats::na.omit(meta$cell_line[match(ids, meta$dataset_id)])))
            }, 0L)
        } else {
            inv <- tryCatch(.mechanism_inventory(drivers, NULL), error = function(e) NULL)
            if (is.null(inv) || !nrow(inv)) integer(0) else
                vapply(drivers, function(g) {
                    d <- inv[inv$perturbed_gene %in% .mechanism_driver_names(g), , drop = FALSE]
                    length(unique(stats::na.omit(d$cell_line)))
                }, 0L)
        }
        ## 用**最好的那个候选**的细胞系数，不用最差的：一个候选覆盖不足只说明
        ## 它自己不该被认领，不代表整份数据不能归属。
        if (length(per)) {
            n_cell_lines <- max(per)
            n_cell_lines_min <- min(per)
        }
    }

    if (!is.na(cnt$n_up) && cnt$n_up < min_up) {
        reasons <- c(reasons, sprintf(
            "few-up-genes: %d genes have log2FC > 1; below %d, pathways are too small to carry a driver claim.",
            cnt$n_up, min_up))
    }
    if (is.finite(med_members) && med_members < min_members) {
        reasons <- c(reasons, sprintf(
            "small-pathways: the median pathway has %d members; below %d there is nothing to test.",
            as.integer(med_members), min_members))
    }
    if (!is.na(n_cell_lines) && n_cell_lines < min_cell_lines) {
        reasons <- c(reasons, sprintf(
            "single-cell-line: the least-covered candidate has data in %d cell line(s); below %d it cannot be separated from an artefact of that line.",
            n_cell_lines, min_cell_lines))
    }
    if (!is.null(reference) && inherits(reference, "mechanismReference")) {
        n_prog <- length(unique(stats::na.omit(as.character(reference$meta$perturbed_gene))))
        if (n_prog < min_background) {
            reasons <- c(reasons, sprintf(
                "thin-reference: only %d comparable programmes; below %d the specificity percentile is unreliable.",
                n_prog, min_background))
        }
    }

    tier <- if (length(reasons)) "description" else "attribution"
    out <- .mechanism_preflight_out(tier, reasons, cnt, med_members, n_cell_lines, cnt$n_up)
    if (exists("n_cell_lines_min")) out$metrics$min_cell_lines <- n_cell_lines_min
    out
}

.mechanism_preflight_out <- function(tier, reasons, cnt, med_members, n_cell_lines, n_up) {
    out <- structure(list(
        tier = tier, reasons = reasons,
        metrics = data.frame(
            n_pathways = cnt$n_pathways, n_up = n_up,
            median_members = med_members, n_cell_lines = n_cell_lines,
            stringsAsFactors = FALSE),
        counts = cnt[c("n_pathways", "n_up")]),
        class = "mechanismPreflight")
    out
}

#' @rdname mechanism_preflight
#' @param ... Ignored.
#' @export
print.mechanismPreflight <- function(x, ...) {
    cat("<mechanism preflight>\n")
    cat("  tier:    ", x$tier, "\n", sep = "")
    m <- x$metrics
    cat(sprintf("  inputs:  %d pathways, %s up genes, median %s members, %s cell line(s)\n",
                m$n_pathways,
                if (is.na(m$n_up)) "?" else as.character(m$n_up),
                if (is.finite(m$median_members)) as.character(as.integer(m$median_members)) else "?",
                if (is.na(m$n_cell_lines)) "?" else as.character(m$n_cell_lines)))
    if (length(x$reasons)) {
        cat("  why:     ", paste(x$reasons, collapse = "\n           "), "\n", sep = "")
    } else {
        cat("  why:      nothing blocked attribution\n")
    }
    invisible(x)
}
