## 通路层面的机制证据。
##
## 与 `evidence()` 的 TF 模式并列，但问的是另一个问题：
##   不是"哪个 TF 调控这批基因"，而是"这条通路的移动能否被一份实测扰动解释"。
##
## 设计依据 `dev_docs/pathway-mechanism-design.md`。两条核心约束：
##
## 1. **候选驱动者由调用方给**（`drivers`）。工具不假装从几千个数据集里发现驱动者——
##    第 0 期已经测过，那种"谁最像"的排序给不出可用答案（biotrigger 把答案埋在第 18/200）。
##    这里做的是**检验**：给一个候选，回答"敲掉它会不会让这条通路动"。
## 2. **每条结论都必须可被推翻**，所以输出带一句 `statement`。不可falsify 的不是机制陈述。

## 通路成员与方向。
##
## 方向优先级：显式的 `pathway_direction` > 从成员基因的 log2FC 推断。
## 两者都没有时方向为 NA——**不猜**。没有方向的通路无法检验，只能给 C 级。
.pathway_members <- function(result) {
    if (is.null(result$geneID)) return(character(0))
    unique(trimws(unlist(strsplit(as.character(result$geneID), "/"))))
}

## 方向从哪来，决定 `support` 能不能当证据用。
##   "given"   —— 结果自己写明（外部信息，support 独立）
##   "nes"     —— GSEA 的 NES 符号（外部信息，support 独立）
##   "members" —— 由成员多数推导（**support 因此是恒真的**，见 .pathway_direction_source）
.pathway_direction_source <- function(result, signed) {
    if (!is.null(result$pathway_direction) && !is.na(result$pathway_direction)) return("given")
    if (!is.null(result$NES) && !is.na(result$NES)) return("nes")
    if (!is.null(signed) && length(signed)) return("members")
    "none"
}

.pathway_direction <- function(result, members, signed) {
    if (!is.null(result$pathway_direction) && !is.na(result$pathway_direction)) {
        return(as.character(result$pathway_direction))
    }
    if (!is.null(result$NES) && !is.na(result$NES)) {
        return(if (result$NES > 0) "up" else "down")
    }
    ## An ORA result carries no direction of its own -- it says the pathway is
    ## enriched, not which way it moved. When the caller supplies signed fold
    ## changes the members answer it: whichever way the majority went is the way
    ## the pathway went. Returning NULL here instead made every row untestable,
    ## which is what the first run of this function did.
    if (is.null(signed) || !length(signed)) return(NULL)
    hit <- members[members %in% names(signed)]
    if (!length(hit)) return(NULL)
    moves <- signed[hit]
    n_up <- sum(moves > 0)
    n_down <- sum(moves < 0)
    if (n_up > n_down) "up" else if (n_down > n_up) "down" else NULL
}

## 成员方向一致性：方向一致且过阈值的成员占比。
##
## 这条存在的理由与旧 spec 的 `support` 相同——少数基因驱动的高分不能算整条通路的支持。
## 一条 33 个成员的通路里若只有 5 个朝同一方向，那不是"通路在动"，是"几个基因在动"。
.pathway_support <- function(members, signed, direction) {
    if (is.null(signed) || !length(signed) || is.null(direction)) return(NA_real_)
    hit <- members[members %in% names(signed)]
    if (!length(hit)) return(NA_real_)
    moves <- signed[hit]
    want <- if (direction == "up") moves > 0 else moves < 0
    mean(want)
}

## 把候选驱动基因解析成 GPSAdb 数据集。
##
## 一个基因可能有多份数据集。全部取回，逐一检验，取最强的那一份——
## 但把用了哪一份记下来，因为"同一个基因在不同细胞系里表现不同"是真实存在的
## （第 0 期里 VHL 在 HK-2 上就不符合预期）。
.pathway_datasets <- function(gene, limit = 5L) {
    d <- tryCatch(as.data.frame(
        gpsadb::gpsadb_search_datasets(perturbed_gene = gene, method = "knockdown",
                                       limit = 50)[[1]]),
        error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- d[!is.na(d$in_reference) & d$in_reference, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    ## 每个细胞系先留一份（样本多的），保证取到的几份尽量来自不同系——
    ## 复现判定看的是细胞系数，不是数据集数。
    d <- d[order(-d$n_samples), , drop = FALSE]
    d <- d[!duplicated(d$cell_line), , drop = FALSE]
    utils::head(d, limit)
}

#' Test enriched pathways against measured perturbation data
#'
#' @param x An `enrichResult` or `gseaResult` from [enrichGO()], [enrichKEGG()]
#'   and friends, or a `$pathways` table from an [evidence()] bundle.
#' @param drivers Character vector of candidate driver genes. Each is resolved to
#'   its perturbation datasets and tested against every pathway.
#' @param signed Optional named numeric of log2 fold changes, used to work out
#'   which way a pathway's members moved when the result does not say.
#' @param alpha Significance threshold for a perturbation to count.
#' @param min_support Minimum fraction of members moving the same way before a
#'   pathway counts as moved at all.
#' @param max_pathways How many pathways to report, taken in the result's own order.
#' @param max_datasets_per_driver Cap on how many datasets are tested per driver.
#'   Testing all of them is `n_drivers * n_datasets` calls and grows without
#'   bound; six is enough to see two cell lines agree, which is all level A
#'   needs. Datasets from distinct cell lines are preferred over extra replicates
#'   of one line, because replication is what the grade turns on.
#' @param n For `print()`, how many rows to show.
#' @param ... Ignored, for method consistency.
#' @return A `data.frame`, one row per pathway, with a `statement` column holding
#'   a claim that a perturbation experiment could refute.
#' @export
pathway_evidence <- function(x, drivers = NULL, signed = NULL, alpha = 0.05,
                             min_support = 0.5, max_pathways = 20L,
                             max_datasets_per_driver = 6L) {
    d <- if (is.data.frame(x)) x else as.data.frame(x)
    if (!nrow(d)) return(.pathway_evidence_empty())
    d <- utils::head(d, max_pathways)

    ## 候选数据集只解析一次，不按通路重复查——同一个驱动者在几十条通路上是同一批数据集。
    sets <- list()
    if (length(drivers)) {
        for (g in drivers) {
            ds <- .pathway_datasets(g, limit = max_datasets_per_driver)
            if (!is.null(ds)) sets[[g]] <- ds[, c("dataset_id", "cell_line"), drop = FALSE]
        }
    }

    rows <- lapply(seq_len(nrow(d)), function(i) {
        members <- .pathway_members(d[i, ])
        dsrc <- .pathway_direction_source(d[i, ], signed)
        direction <- .pathway_direction(d[i, ], members, signed)
        support <- .pathway_support(members, signed, direction)
        desc <- d$Description[i] %||% d$ID[i] %||% NA_character_
        ## 先给默认值：下面的"无方向"分支会提前 return，那时检验还没跑，
        ## 但已经要填 `n_datasets_tested` 这一列。
        n_tested <- 0L

        ## 没有方向有两种，必须是两个结论，不是一个：
        ##   * `signed` 没给 —— 我们没有测过成员往哪边走，这是**没测**；
        ##   * `signed` 给了但成员对半分 —— 测过了，它们没有共同方向，这是**测了不成立**。
        ## 混成一句"没有方向"会把后者说成缺失，与 `absent` / `untested` 分开是同一条纪律。
        if (is.null(direction)) {
            untested <- is.null(signed) || !length(signed)
            return(data.frame(
                pathway = desc, direction = NA_character_, n_members = length(members),
                support = support, driver = NA_character_, dataset = NA_character_,
                cell_line = NA_character_, es = NA_real_, padj = NA_real_,
                effect = NA_character_, n_datasets_tested = n_tested,
                n_cell_lines_supported = 0L, level = "C",
                verdict = if (untested) "untestable" else "unsupported",
                statement = if (untested)
                    "No signed fold changes were supplied, so the pathway's movement is unknown."
                else
                    "Members split evenly, so the pathway has no common direction to explain.",
                stringsAsFactors = FALSE))
        }

        ## 逐个候选驱动者的**每一份**数据集都检验，不只看最强的那一份。
        ##
        ## 原版只留 padj 最小的一个，于是 A 级等于"某一份数据支持"。实测这个含义太弱：
        ## 四个驱动者里三个换一个细胞系就不成立，而且效应量符号翻转
        ## （`results/stage3_refutation.tsv`）。所以要看得见全部数据集，才能判断复现。
        tested <- list()
        for (g in names(sets)) {
            for (j in seq_len(nrow(sets[[g]]))) {
                id <- as.character(sets[[g]]$dataset_id[j])
                r <- tryCatch(gpsadb::gpsadb_test_gene_set(id, members),
                              error = function(e) NULL)
                if (is.null(r)) next
                t <- as.data.frame(r[[1]])
                if (!nrow(t) || is.na(t$padj[1])) next
                tested[[length(tested) + 1L]] <- list(
                    gene = g, dataset = id,
                    cell_line = as.character(sets[[g]]$cell_line[j]),
                    es = as.numeric(t$es[1]), padj = as.numeric(t$padj[1]))
            }
        }
        ## 最强的显著者决定"提出主张"的那一份；复现则要求**另一个细胞系**也显著且同向。
        best <- NULL
        for (cand in tested) {
            if (is.na(cand$padj) || cand$padj >= alpha) next
            if (is.null(best) || cand$padj < best$padj) best <- cand
        }
        n_tested <- length(tested)
        cells_agree <- character(0)
        if (!is.null(best)) {
            want <- sign(best$es)
            cells_agree <- unique(vapply(tested, function(cand) {
                if (is.na(cand$padj) || cand$padj >= alpha) return(NA_character_)
                if (is.na(cand$es) || sign(cand$es) != want) return(NA_character_)
                cand$cell_line
            }, character(1)))
            cells_agree <- cells_agree[!is.na(cells_agree)]
        }
        replicated <- length(cells_agree) >= 2L

        if (!is.null(best) && best$padj < alpha) {
            ## 显著性还不够：**方向必须对得上**。
            ## 通路往上走、驱动者被敲低，那这条通路应当**回落**（es < 0）——这说明
            ## 该因子在启动它。若敲低后通路反而**更往上**（es > 0），那它是在压制，
            ## 不是驱动。原版不看这一点，把任何显著结果都记成支持。
            effect <- if (is.na(best$es) || is.na(direction)) NA_character_ else
                if ((direction == "up" && best$es < 0) ||
                    (direction == "down" && best$es > 0)) "activator" else "repressor"
            agrees <- is.na(effect) || effect == "activator"
            ## A 级要求**复现**：至少两个不同细胞系都显著且同向。
            ## 单份数据支持降为 B 并标记 `single-dataset`——它不是没有证据，
            ## 是证据只有一个来源，而实测这种证据有 3/4 会在换数据后翻转。
            level <- if (!agrees) "B" else if (replicated) "A" else "B"
            verdict <- if (!agrees) "direction-conflict"
                       else if (replicated) "perturbation-supported"
                       else "single-dataset"
            statement <- if (agrees && replicated) sprintf(
                "If %s drives this programme, knocking it down should move these %d members; %s does (padj %.2g) and %d cell lines agree.",
                best$gene, length(members), best$dataset, best$padj, length(cells_agree))
            else if (agrees) sprintf(
                "If %s drives this programme, knocking it down should move these %d members; %s does (padj %.2g), but only one cell line supports it.",
                best$gene, length(members), best$dataset, best$padj)
            else sprintf(
                "%s moves these %d members in the same direction the pathway already runs (padj %.2g), which reads as repression rather than drive.",
                best$gene, length(members), best$padj)
            return(data.frame(
                pathway = desc, direction = direction, n_members = length(members),
                support = support, driver = best$gene, dataset = best$dataset,
                cell_line = best$cell_line, es = best$es, padj = best$padj,
                effect = effect, n_datasets_tested = n_tested,
                n_cell_lines_supported = length(cells_agree),
                level = level, verdict = verdict, statement = statement,
                stringsAsFactors = FALSE))
        }

        ## 没有实测支持时，成员方向一致只能算弱证据（B 级）——而且**只有当方向来自外部时**才算。
        ##
        ## 方向若由成员多数推导，`support` 就是恒真的（多数必然 >= 50%），拿它当证据
        ## 等于用同一批基因既定方向又证方向。这是审核时发现的：`min_support = 0.5`
        ## 永远成立，B 级因此没有信息量。
        if (!is.na(support) && support >= min_support && dsrc != "members") {
            return(data.frame(
                pathway = desc, direction = direction, n_members = length(members),
                support = support, driver = NA_character_, dataset = NA_character_,
                cell_line = NA_character_, es = NA_real_, padj = NA_real_,
                effect = NA_character_, n_datasets_tested = n_tested,
                n_cell_lines_supported = 0L, level = "B", verdict = "coherent-untested",
                statement = sprintf(
                    "%.0f%% of members moved %s, so the pathway moved; no perturbation test identified a driver.",
                    100 * support, direction),
                stringsAsFactors = FALSE))
        }

        data.frame(
            pathway = desc, direction = direction, n_members = length(members),
            support = support, driver = NA_character_, dataset = NA_character_,
            cell_line = NA_character_, es = NA_real_, padj = NA_real_,
            effect = NA_character_, n_datasets_tested = n_tested,
            n_cell_lines_supported = 0L, level = "C", verdict = "unsupported",
            statement = if (is.na(support))
                "Members were not found in the signed gene list; the pathway's movement is unknown."
            else sprintf("%.0f%% of members moved %s, below the support floor of %.0f%%.",
                         100 * support, direction, 100 * min_support),
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, rows)
    out[order(match(out$level, c("A", "B", "C")), -out$support), , drop = FALSE]
}

.pathway_evidence_empty <- function() {
    data.frame(pathway = character(), direction = character(), n_members = integer(),
               support = numeric(), driver = character(), dataset = character(),
               cell_line = character(), es = numeric(), padj = numeric(),
               effect = character(), n_datasets_tested = integer(),
               n_cell_lines_supported = integer(),
               level = character(), verdict = character(), statement = character(),
               stringsAsFactors = FALSE)
}

#' @rdname pathway_evidence
#' @export
print.pathway_evidence <- function(x, n = 10L, ...) {
    cat("<pathway evidence>\n")
    counts <- table(factor(x$level, levels = c("A", "B", "C")))
    cat("  levels:       ", paste0(names(counts), " ", as.integer(counts),
                                  collapse = "  "), "\n", sep = "")

    shown <- utils::head(x, n)
    if (nrow(shown)) {
        cat("\n")
        for (i in seq_len(nrow(shown))) {
            row <- shown[i, ]
            driver <- if (is.na(row$driver)) "-" else
                sprintf("%s/%s padj %.2g", row$driver, row$cell_line, row$padj)
            cat(sprintf("  [%s] %s (%s, %d members, %.0f%% together)\n",
                        row$level, substr(row$pathway, 1, 52), row$direction %||% "?",
                        row$n_members, 100 * (row$support %||% NA)))
            cat(sprintf("      driver: %s\n", driver))
        }
        if (nrow(x) > nrow(shown)) {
            cat(sprintf("  ... %d more\n", nrow(x) - nrow(shown)))
        }
    }

    ## 陈述放在末尾而不是每行都印：它是要被检验的主张，不是摘要的一部分。
    a <- x[x$level == "A", , drop = FALSE]
    if (nrow(a)) {
        cat("\n  claim to test:\n")
        cat("   ", a$statement[1], "\n")
    }
    invisible(x)
}
