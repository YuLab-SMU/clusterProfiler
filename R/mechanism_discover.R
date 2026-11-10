## 发现门：从一条差异表达谱，排出"可能驱动它的转录因子"。
##
## 为什么要有这一层。核验层（`pathway_evidence()`）回答的是"给定一个候选，证据够不够"，
## 它按设计不能回答"是谁"——每份数据要 265 个候选全测一遍才能识别，按每次 GPSAdb
## 调用 1–2 秒算是 30–50 分钟一份。所以要有一个便宜的第一道：先排出候选，再逐个核验。
##
## 方法不是我们发明的，是 `ulm`（decoupler）：把某个 TF 的靶基因权重向量（靶 = ±1，
## 非靶 = 0）与整条 log2FC 谱求 Pearson 相关，再用样本量换成 t 值。
## 我们试过 12 个组合（集成、跨网络、viper、按绝对值排序），没有一个超过它。
##
## Full-universe benchmark and frozen inputs: dev_docs/mechanism-benchmark-4215.R.
## A missing top-20 rank is a ranking miss, not evidence of network absence.
##
## 所以这里不"改进"它，只做三件它没有的事：**去自环、给方向、说清覆盖多少**。

#' Rank candidate transcription factors for a differential expression profile
#'
#' This is the cheap first stage. It ranks transcription factors by correlating
#' each one's signed target weights against the whole differential expression
#' vector, then leaves the decision to [pathway_evidence()], which checks a
#' named candidate against measured perturbation data.
#'
#' The score is the t-value of that correlation. **A positive score means the
#' factor looks active in this experiment, a negative one means repressed.**
#' The ranking is not a verdict. In a frozen exploratory corpus of 17 profiles,
#' all 17 expected drivers were present and scorable with CollecTRI. Using the
#' default signed ranking, expected-label hit counts were 6/17 at rank one,
#' 10/17 in the top five and 13/17 in the top ten. The corpus includes one
#' knockout and repeated studies/drivers; these are not independent validation
#' estimates or proof of unique causal identification. A rank outside the first
#' twenty must not be counted as network absence. Treat the output as a shortlist.
#'
#' In the truncated top-twenty candidate pool, score and measured target count
#' were positively associated, including among non-reference candidates
#' (Spearman correlation about 0.24). This selected-pool observation does not
#' establish a general size bias or show whether annotation depth, correlated
#' targets or biology explains it. `n_targets` and `share` disclose coverage.
#' No target-count weight is applied: exploratory adjustments have not been
#' independently validated.
#'
#' @param signed Named numeric vector of log2 fold changes, names are gene
#'   symbols. The whole vector is used, not just the significant genes: the
#'   correlation is what separates a factor's targets from everything else.
#' @param network Data frame of transcription factor targets with columns
#'   `source` (factor), `target` (gene symbol) and optionally `weight`
#'   (positive for activation targets, negative for repression targets; 1 when
#'   absent). CollecTRI and DoRothEA can be exported in this shape. Identical
#'   repeated edges are deduplicated; conflicting weights for the same
#'   source-target pair are rejected rather than resolved by row order.
#' @param top Integer. How many candidates to return.
#' @param min_targets Integer. Factors with fewer measured targets than this are
#'   skipped: a correlation over three genes is noise.
#' @param rank_by Character. `"score"` (default) ranks by the signed score, so
#'   factors that look active come first; that is the right default for
#'   treatment and ligand experiments, where the driver is switched on.
#'   `"abs"` ranks by magnitude regardless of sign, which is what a knockdown
#'   or knockout profile needs: its driver scores negative and would otherwise
#'   sink to the bottom. On our benchmark 16 of 17 profiles are activation
#'   types, so `"score"` wins there; that is a property of the corpus, not of
#'   the method.
#' @param keep_self_targets Logical. When `FALSE` (the default) a factor that
#'   lists itself as a target loses that edge. Self-targets inflate the score
#'   for any knockdown or knockout, because the factor's own gene is the
#'   strongest signal in the profile; in CollecTRI 62% of the factors that
#'   reach evaluation carry such an edge, and removing them costs about five
#'   percentage points of ranking performance.
#'
#' @return A data frame with one row per candidate: `rank`, `source`, `score`
#'   (the t-value), `r` (the correlation), `direction` (`"activates"` or
#'   `"represses"`), `n_targets`, and `share` (the fraction of all measured
#'   genes that are targets of this factor; the denominator is `length(signed)`,
#'   not the number of genes in the network).
#'
#' @seealso [pathway_evidence()] to test the shortlist against perturbation
#'   data, and [mechanism_preflight()] to decide whether the profile can carry
#'   any claim at all.
#'
#' @examples
#' set.seed(2)
#' genes <- paste0("g", 1:50)
#' signed <- stats::setNames(c(rep(1.5, 5), stats::rnorm(45, 0, 0.2)), genes)
#' network <- data.frame(source = "TF1", target = genes[1:5], weight = 1)
#' mechanism_discover(signed, network)
#'
#' @export
mechanism_discover <- function(signed, network, top = 10L, min_targets = 5L,
                               rank_by = c("score", "abs"),
                               keep_self_targets = FALSE) {
    rank_by <- match.arg(rank_by)
    if (!is.numeric(signed) || is.null(names(signed)))
        stop("`signed` must be a named numeric vector (gene symbol to log2 fold change).",
             call. = FALSE)
    ## 重复基因名会让 `w[match()]`（只填第一个位置）与 `x`（每个重复位置都有值）
    ## 错位，相关分数静默给错。同包的 `finish()` 早就按 |lfc| 最大去重，这里照做。
    if (anyDuplicated(names(signed))) {
        n_dup <- sum(duplicated(names(signed)))
        ord <- order(-abs(signed))
        signed <- signed[ord]
        signed <- signed[!duplicated(names(signed))]
        message(sprintf("Dropped %d duplicated gene name(s), keeping the largest |log2FC|.", n_dup))
    }
    net <- as.data.frame(network, stringsAsFactors = FALSE)
    if (!all(c("source", "target") %in% names(net)))
        stop("`network` must have `source` and `target` columns.", call. = FALSE)
    if (!"weight" %in% names(net)) net$weight <- 1
    net$source <- as.character(net$source)
    net$target <- as.character(net$target)
    net$weight <- as.numeric(as.character(net$weight))
    net <- net[!is.na(net$source) & !is.na(net$target) & is.finite(net$weight), ,
               drop = FALSE]

    ## 自环：靶基因里含 TF 自己。不删的话，敲低/敲除实验里 TF 自身的 log2FC
    ## 会直接把它的分数顶上去（CollecTRI 里 62% 的 TF 带这种边）。
    if (!keep_self_targets) {
        self <- net$source == net$target
        if (any(self)) {
            message(sprintf("Dropped %d self-target edges (%d factor(s)). See `keep_self_targets`.",
                            sum(self), length(unique(net$source[self]))))
        }
        net <- net[!self, , drop = FALSE]
    }
    if (!nrow(net))
        stop("No usable edges left in `network`.", call. = FALSE)

    ## NA 的处理：`signed` 里可能有没算出来的基因。留着会让整列相关变成 NA，
    ## 而"这条基因没测出 fold change"与"它没变化"是两件事，所以直接剔除并说明。
    if (any(!is.finite(signed))) {
        n_na <- sum(!is.finite(signed))
        message(sprintf("Dropped %d gene(s) with non-finite log2 fold change.", n_na))
        signed <- signed[is.finite(signed)]
    }
    if (!length(signed)) stop("`signed` holds no finite values.", call. = FALSE)

    ## Resolve exact duplicates only; conflicting evidence must not depend on row order.
    net <- net[!duplicated(net[, c("source", "target", "weight")]), , drop = FALSE]
    if (anyDuplicated(net[, c("source", "target")]))
        stop("Conflicting weights for the same source-target pair in `network`; resolve them before scoring.",
             call. = FALSE)
    net <- net[net$target %in% names(signed), , drop = FALSE]
    if (!nrow(net))
        stop("None of the network targets are measured in `signed`.", call. = FALSE)

    ## 所有 TF 用**同一批基因**算相关，df 才是常数。
    ##
    ## 基因域取**整条测量谱**，不是只取网络靶基因。两个理由：
    ##   1. 统计上要的是"这组靶基因相对全转录组特殊吗"，只取靶基因就没有对照组了
    ##      （单 TF 网络里 w 全是 1，方差为 0，函数会直接拒绝）；
    ##   2. 网络总是残缺的，没被任何 TF 收录的基因仍是背景的一部分。
    ## 网络覆盖广时（CollecTRI 约 4600 个靶），这与 decoupler 在 adj 特征上算的结果一致。
    genes <- names(signed)
    x <- as.numeric(signed[genes])
    xc <- x - mean(x)
    sx <- sqrt(sum(xc^2))
    if (!is.finite(sx) || sx == 0)
        stop("The measured part of `signed` has no variance.", call. = FALSE)

    split_idx <- split(seq_len(nrow(net)), net$source)
    rows <- lapply(names(split_idx), function(tf) {
        idx <- split_idx[[tf]]
        w <- numeric(length(genes))
        w[match(net$target[idx], genes)] <- net$weight[idx]
        n <- sum(w != 0)
        if (n < min_targets) return(NULL)
        wc <- w - mean(w)
        sw <- sqrt(sum(wc^2))
        if (!is.finite(sw) || sw == 0) return(NULL)
        r <- sum(wc * xc) / (sw * sx)
        r <- max(min(r, 1 - 1e-12), -1 + 1e-12)
        df <- length(genes) - 2L
        data.frame(source = tf, r = r,
                   score = r * sqrt(df / ((1 - r) * (1 + r))),
                   direction = if (r > 0) "activates" else "represses",
                   n_targets = n, share = n / length(genes),
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    if (is.null(out) || !nrow(out))
        stop("No factor had enough measured targets. Lower `min_targets` or widen the network.",
             call. = FALSE)
    ## 实际参与评分的因子数，不是网络里有边的因子数。后者会把因 `min_targets`
    ## 被跳过的也算进去，用户看到 "tested 1176" 会以为全都评了。
    n_scored <- sum(!vapply(rows, is.null, logical(1)))
    out <- out[order(if (rank_by == "abs") -abs(out$score) else -out$score), ,
               drop = FALSE]
    out$rank <- seq_len(nrow(out))
    out <- utils::head(out, max(1L, as.integer(top)[1]))
    out <- out[, c("rank", "source", "score", "r", "direction", "n_targets",
                   "share"), drop = FALSE]
    rownames(out) <- NULL
    ## 属性在取子集之后设置：`[.data.frame` 会丢掉它们，放在前面就白设了。
    attr(out, "n_genes") <- length(genes)
    attr(out, "n_factors_tested") <- n_scored
    attr(out, "n_factors_in_network") <- length(split_idx)
    class(out) <- c("mechanism_discovery", "data.frame")
    out
}

#' @export
print.mechanism_discovery <- function(x, n = 10L, ...) {
    cat("<mechanism discovery>\n")
    cat("  genes used:  ", attr(x, "n_genes"), "\n", sep = "")
    cat("  factors tested:", attr(x, "n_factors_tested"),
        " (of ", attr(x, "n_factors_in_network"), " in network)\n", sep = "")
    show <- utils::head(as.data.frame(x), n)
    print(show[, c("rank", "source", "score", "direction", "n_targets")],
          row.names = FALSE)
    cat("\nThis is a shortlist, not a verdict: verify the candidates with pathway_evidence().\n")
    invisible(x)
}
