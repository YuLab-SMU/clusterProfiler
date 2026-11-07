## 判决书层：把「结论」换成「已排除什么 / 还缺什么 / 下一次做什么」。
##
## 为什么要有这一层。九次实跑暴露的事实是不对称的：这套方法**排除**很准
## （12 个不相干的驱动者里 5 个被判 `direction-conflict`，全对），**识别**很弱
## （驱动者层面 13% 假阳性，全是必需基因）。在这种不对称下，"每行给一个等级"
## 是最不容易被读对的形式——读者看到一排 C，分不清"没找到"与"没法找"，
## 而实测十一个数据集里七个是后者。
##
## 所以把 `reason` 从"降级说明"提升为**主产品**：每条通路一行判决，
## 每个缺口一条可执行的下一次实验，并且**允许"没有实验能救"这个答案**——
## 差异基因太散、通路太小这类问题，再补一份扰动数据也不会变。

## 补救动作表。
##
## 三侧：data（补数据）、analysis（改输入或口径）、concept（这份数据不适合这套方法）。
## 第三侧必须存在。把 analysis 与 concept 的失败写成"再补一份扰动数据"，
## 是让用户花钱去做一件不会改变结论的事。
.mechanism_remedies <- function() {
    ## 一行一条，四个字段写在一起。原来写四个平行向量，加一条 `no-background`
    ## 时 `side` 少补一个，`data.frame()` 当场报 "differing number of rows: 18, 17"
    ## ——位置对齐的错误在增删条目时必然复发，所以换成按行声明。
    rows <- list(
        list("single-cell-line", "data", "add a knockdown dataset in a second cell line", TRUE),
        list("no-dataset", "data", "no perturbation data exists for this driver; a knockdown or overexpression in any cell line would start it", TRUE),
        list("no-test-completed", "data", "the calls failed; rerun before reading these rows", TRUE),
        list("too-small", "analysis", "supply more differentially expressed genes, or use a library whose pathways are larger", TRUE),
        list("low-support", "analysis", "members moved inconsistently; widen the contrast or lower the |log2FC| floor", FALSE),
        list("member-derived", "analysis", "supply fold changes so the direction is measured rather than inferred from the same members", FALSE),
        list("no-members", "analysis", "the pathway's members are missing from the signed list; check the identifier space", FALSE),
        list("no-fold-changes", "analysis", "no fold changes were supplied", FALSE),
        list("no-direction", "analysis", "the members have no common direction", FALSE),
        list("not-significant", "concept", "nothing to add: the test ran and the driver did not move this pathway", FALSE),
        list("coherent-refuted", "concept", "nothing to add: the tests ran and none moved this pathway", FALSE),
        list("no-perturbation-test", "analysis", "no test could run without a driver and its data", FALSE),
        ## direction-conflict 是**证据反对**（敲低后通路同向移动，读作压制），
        ## 不是"数据不够"。归 data 侧且 unlocks=TRUE 会承诺"再补一份数据就能救"，
        ## 与判决书允许"没有实验能救"的设计原则相悖。
        list("direction-conflict", "concept", "nothing to add here: the perturbation moved this pathway with the condition, so the candidate reads as repression rather than drive", FALSE),
        list("supported", "concept", "validate the driver in an independent perturbation; an A row is consistency worth following up, not proof", FALSE))
    data.frame(
        reason = vapply(rows, function(z) z[[1]], character(1)),
        side = vapply(rows, function(z) z[[2]], character(1)),
        action = vapply(rows, function(z) z[[3]], character(1)),
        unlocks = vapply(rows, function(z) isTRUE(z[[4]]), logical(1)),
        stringsAsFactors = FALSE)
}

.mechanism_remedy_for <- function(reason, verdict) {
    tab <- .mechanism_remedies()
    key <- reason
    if (is.na(key) || !nzchar(key)) {
        ## 等级 A 与 `specific` 是同一种"已建立"，下一步也应当是同一句：
        ## 在独立扰动里验证。只给其中一个写动作，读者会以为另一个不需要验证。
        ## `verdict == "too-small"` 在 pathway_evidence 里目前不可达，但键必须是
        ## 补救表里有的那一个：原来写 `"too-small-decomp"`（全仓只此一处），
        ## match 落空后静默返回 action = NA，那一行就没有下一步了。
        key <- if (identical(verdict, "perturbation-supported")) "supported"
            else if (identical(verdict, "too-small")) "too-small"
            else NA_character_
    }
    hit <- tab[match(key, tab$reason), , drop = FALSE]
    if (!nrow(hit) || is.na(hit$side)) {
        return(list(side = NA_character_, action = NA_character_, unlocks = FALSE))
    }
    list(side = hit$side, action = hit$action, unlocks = hit$unlocks)
}

## 每条通路落在三档里的哪一档。
##
## 措辞从 established/excluded/undecided 改成 supported/excluded/undetermined：
## 前一组读起来像"结论已经建立"，那是这条通路给不了的；后一组只承诺
## "有 A 级证据 / 有反证或测过不支持 / 还看不清"。
##
## supported   —— 有 A 级证据（至少两个细胞系同向）；
## excluded    —— 看了，没有程序解释它，或证据方向相反；
## undetermined —— 通路太小、没有覆盖，或缺输入。**这一档才是报告的主体。**
.mechanism_verdict_class <- function(x) {
    if ("level" %in% names(x)) {
        return(ifelse(x$level == "A", "supported",
               ifelse(x$verdict %in% c("direction-conflict", "coherent-refuted"), "excluded",
                      "undetermined")))
    }
    rep("undetermined", nrow(x))
}

## 缺口：该驱动者在库里**全部**可用数据集，而不只是被抽样测过的那几份。
##
## 这一点是产品缺口：现在 bundle 里只有抽到的六份，于是 `reason` 说得出"不够"，
## 说不出"还差什么"。要给出可执行的下一步，必须看见完整清单。
## 别名候选。与 `pathway_evidence()` 那边 `.pathway_datasets()` 用的是同一张表、
## 同一个双向查找：GPSAdb 里 `GR` 没有数据、`NR3C1` 有，而用户写 `drivers = "GR"`
## 是最自然的。不做这一步，报告会对 `GR` 说"库里没有这个驱动者的任何数据"——
## 一句**假的**事实陈述，而不是"没查到"。
.mechanism_driver_names <- function(gene) {
    aliases <- tryCatch(mechanism_regulator_aliases(), error = function(e) character(0))
    back <- if (length(aliases)) names(aliases)[aliases == gene] else character(0)
    fwd <- tryCatch(mechanism_alias_regulator(gene), error = function(e) gene)
    tries <- unique(c(gene, fwd, back))
    tries[!is.na(tries) & nzchar(tries)]
}

.mechanism_inventory <- function(drivers, inventory = NULL) {
    if (is.null(inventory)) {
        if (!length(drivers)) return(NULL)
        if (!requireNamespace("gpsadb", quietly = TRUE)) return(NULL)
        parts <- lapply(drivers, function(g) {
            found <- NULL
            for (nm in .mechanism_driver_names(g)) {
                d <- tryCatch(as.data.frame(
                    gpsadb::gpsadb_search_datasets(perturbed_gene = nm, method = "knockdown",
                                                   limit = 100)[[1]]),
                    error = function(e) NULL)
                if (!is.null(d) && nrow(d)) { found <- d; found$queried_as <- g; break }
            }
            found
        })
        inventory <- do.call(rbind, parts[!vapply(parts, is.null, TRUE)])
        if (is.null(inventory) || !nrow(inventory)) return(NULL)
    }
    inventory <- as.data.frame(inventory, stringsAsFactors = FALSE)
    ## 缺 `perturbed_gene` 时不能当成"这个驱动者没有数据"——那是把调用方传错列名
    ## 说成一条生物学结论。实测把一个叫 `driver` 的列传进来，四个驱动者全被报成
    ## "no dataset in GPSAdb"，而库里明明各有 2-6 份。
    if (!"perturbed_gene" %in% names(inventory)) {
        rlang::abort(c("`inventory` has no `perturbed_gene` column.",
                       i = "That column says which gene each dataset knocked down; without it every driver reads as having no data.",
                       i = "Rename the column, or leave `inventory = NULL` to query GPSAdb."))
    }
    keep <- intersect(c("dataset_id", "perturbed_gene", "cell_line", "tissue",
                        "method", "n_samples", "in_reference", "queried_as"),
                      names(inventory))
    inventory <- inventory[, keep, drop = FALSE]
    ## 同一份数据集只能算一次。别名查询（`GR` 与 `NR3C1` 指向同一批数据）或者调用方
    ## 把两次查询的结果拼起来时，同一行会出现两遍，于是"1 份数据、1 个细胞系"被数成
    ## 2 份——而"够不够两个细胞系"正是这个领域里的一条门槛。
    if ("dataset_id" %in% names(inventory)) {
        inventory <- inventory[!duplicated(inventory$dataset_id), , drop = FALSE]
    }
    inventory
}

.mechanism_gap_table <- function(inventory, drivers, context) {
    if (is.null(inventory) || !nrow(inventory)) return(NULL)
    tissue <- if (is.list(context)) context$tissue else NULL
    rows <- lapply(drivers, function(g) {
        ## 查询名与别名都算，`GR` 与 `NR3C1` 是同一个驱动者的两种写法。
        d <- inventory[inventory$perturbed_gene %in% .mechanism_driver_names(g), ,
                       drop = FALSE]
        if (!nrow(d)) {
            return(data.frame(driver = g, n_datasets = 0L, n_cell_lines = 0L,
                              tissues = NA_character_, in_context_tissue = 0L,
                              gap = "no dataset in GPSAdb", stringsAsFactors = FALSE))
        }
        cells <- if ("cell_line" %in% names(d)) unique(stats::na.omit(d$cell_line)) else character(0)
        tissues <- if ("tissue" %in% names(d)) unique(stats::na.omit(d$tissue)) else character(0)
        matched <- if (!is.null(tissue) && length(tissues))
            sum(tolower(tissues) == tolower(tissue)) else NA_integer_
        gap <- if (!is.null(tissue) && !is.na(matched) && matched == 0L)
            sprintf("no dataset from %s; available tissues: %s", tissue,
                    paste(utils::head(tissues, 5), collapse = ", "))
        else if (length(cells) < 2L)
            sprintf("one cell line only (%s)", paste(cells, collapse = ", "))
        else NA_character_
        data.frame(driver = g, n_datasets = nrow(d), n_cell_lines = length(cells),
                   tissues = paste(utils::head(tissues, 5), collapse = ", "),
                   in_context_tissue = matched, gap = gap, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' Turn mechanism evidence into a verdict sheet and a next experiment
#'
#' Report what has been established, what has been ruled out, and what is simply
#' unknown -- and, for each gap, the cheapest thing that would close it.
#'
#' The asymmetry this is built around is measured, not assumed: on the hypoxia
#' pilot the perturbation test ruled out the right candidates five times out of
#' five and admitted an unrelated essential gene thirteen per cent of the time.
#' A tool with that profile is a good exclusion instrument. Reporting a row per
#' pathway with a letter grade hides it, because a table of `C` reads the same
#' whether the test ran and failed or could never run at all -- and in eleven
#' real datasets seven were the second case.
#'
#' @param x A table from [pathway_evidence()].
#' @param drivers Drivers to inventory. Defaults to the drivers named in `x`.
#' @param inventory Optional data.frame of available datasets, with columns like
#'   `perturbed_gene`, `cell_line`, `tissue`. When omitted the whole inventory is
#'   read from GPSAdb, which is the point: the tested rows are a sample of six
#'   datasets per driver, so they can say "not enough" but never "what is
#'   missing".
#' @param preflight Optional result of [mechanism_preflight()]. When it says
#'   attribution is not supported, the report stops recommending experiments:
#'   the honest answer is that no experiment would change these rows, and saying
#'   "add another cell line" to a dataset with 23 genes is advice that cannot pay
#'   off.
#' @param context Named list describing the user's experiment. `context$tissue`
#'   is compared against the tissues the available datasets were run in, so the
#'   report can say "no dataset from airway" instead of "add another cell line".
#' @return A `mechanismReport`: a list with `verdicts` (one row per pathway with
#'   `class` in established/excluded/undecided and the action that follows),
#'   `gaps` (one row per driver), `next_steps` (actions ranked by how many rows
#'   they would unlock) and `summary` (two sentences).
#' @seealso [pathway_evidence()], [mechanism_preflight()]
#' @export
mechanism_report <- function(x, drivers = NULL, inventory = NULL,
                             context = list(), preflight = NULL) {
    if (!is.data.frame(x) || !nrow(x)) {
        return(.mechanism_report_empty(
            "no pathway was scored, so there is nothing to judge; supply more differentially expressed genes or a library with larger pathways."))
    }
    if (!"level" %in% names(x)) {
        rlang::abort(c("`x` has no `level` column, so its rows cannot be classified.",
                       i = "Pass the full table from pathway_evidence(); a column subset drops what this report reads."))
    }
    if (is.null(drivers)) {
        drivers <- unique(stats::na.omit(x$driver))
        if (!length(drivers)) drivers <- attr(x, "drivers")
        if (!length(drivers)) drivers <- character(0)
    }
    drivers <- unique(as.character(drivers))

    klass <- .mechanism_verdict_class(x)
    reason <- if ("reason" %in% names(x)) x$reason else NA_character_
    ## 描述层没有 `verdict` 列，用映射出来的档位顶上——`NULL[i]` 会让后面的
    ## `%in%` 得到 logical(0)，`if` 直接报 "argument is of length zero"。
    vcol <- if ("verdict" %in% names(x)) x$verdict else klass
    rem <- lapply(seq_len(nrow(x)), function(i) .mechanism_remedy_for(reason[i], vcol[i]))

    verdicts <- data.frame(
        pathway = x$pathway,
        class = klass,
        verdict = if ("verdict" %in% names(x)) x$verdict else klass,
        reason = reason,
        detail = if ("driver" %in% names(x)) x$driver else NA_character_,
        side = vapply(rem, function(z) z$side %||% NA_character_, character(1)),
        action = vapply(rem, function(z) z$action %||% NA_character_, character(1)),
        stringsAsFactors = FALSE)

    inv <- .mechanism_inventory(drivers, inventory)
    gaps <- .mechanism_gap_table(inv, drivers, context)

    ## 下一步按"能解锁多少行"排序。只排 data / analysis 两侧——concept 那一侧的
    ## 动作是"没有动作"，混进排序会让它看起来像可执行项。
    ## 预检门说"不支持归属"时，不再推荐实验。
    ## 一份 23 个差异基因的数据，补多少个细胞系都不会让它能归属——那类建议
    ## 只会让用户花钱去做一件不会改变结论的事。
    blocked <- !is.null(preflight) && inherits(preflight, "mechanismPreflight") &&
        !identical(preflight$tier, "attribution")
    actionable <- if (blocked) verdicts[0, , drop = FALSE] else
        verdicts[!is.na(verdicts$action) &
                 !is.na(verdicts$side) & verdicts$side != "concept", , drop = FALSE]
    next_steps <- if (!nrow(actionable)) {
        data.frame(action = character(), side = character(), unlocks = integer(),
                   examples = character(), stringsAsFactors = FALSE)
    } else {
        agg <- lapply(split(actionable, actionable$action), function(z)
            data.frame(action = z$action[1], side = z$side[1], unlocks = nrow(z),
                       examples = paste(utils::head(z$pathway, 3), collapse = "; "),
                       stringsAsFactors = FALSE))
        agg <- do.call(rbind, agg)
        agg[order(-agg$unlocks), , drop = FALSE]
    }

    n <- nrow(verdicts)
    n_est <- sum(verdicts$class == "supported")
    n_exc <- sum(verdicts$class == "excluded")
    n_und <- sum(verdicts$class == "undetermined")
    worst <- if (n_und && !is.null(gaps) && nrow(gaps)) gaps$gap[!is.na(gaps$gap)][1] else NA_character_
    summary <- c(
        if (blocked) sprintf("attribution is not supported by this data (%s); the rows below describe what moved, not who drove it.",
                             paste(sub(":.*", "", preflight$reasons), collapse = ", ")),
        sprintf("%d pathways: %d supported, %d excluded, %d undetermined.", n, n_est, n_exc, n_und),
        if (!is.na(worst)) paste0("Largest gap: ", worst, ".") else NA_character_,
        if (blocked) "No experiment would change these rows: the limit is the input, not the evidence."
        else if (nrow(next_steps)) paste0("Next: ", next_steps$action[1], " (unlocks ",
                                     next_steps$unlocks[1], " row(s)).") else
            "No experiment would change these verdicts.")

    out <- structure(list(
        verdicts = verdicts, gaps = gaps, next_steps = next_steps,
        summary = summary[!is.na(summary)],
        context = context, inventory = inv),
        class = "mechanismReport")
    out
}

.mechanism_report_empty <- function(reason) {
    out <- structure(list(
        verdicts = data.frame(pathway = character(), class = character(),
                              verdict = character(), reason = character(),
                              detail = character(), side = character(),
                              action = character(), stringsAsFactors = FALSE),
        gaps = NULL,
        next_steps = data.frame(action = character(), side = character(),
                                unlocks = integer(), examples = character(),
                                stringsAsFactors = FALSE),
        summary = c(reason, "No experiment would change this: the limit is the input."),
        context = list(), inventory = NULL),
        class = "mechanismReport")
    out
}

#' @rdname mechanism_report
#' @param n Number of rows to show.
#' @param ... Ignored.
#' @export
print.mechanismReport <- function(x, n = 6L, ...) {
    cat("<mechanism report>\n")
    for (s in x$summary) cat("  ", s, "\n", sep = "")
    v <- x$verdicts
    if (nrow(v)) {
        cat("\n  supported:\n")
        e <- utils::head(v[v$class == "supported", , drop = FALSE], n)
        if (nrow(e)) {
            for (i in seq_len(nrow(e))) cat(sprintf("    - %s <- %s\n", substr(e$pathway[i], 1, 56), e$detail[i]))
        } else cat("    (none)\n")
        cat("\n  excluded:\n")
        ex <- utils::head(v[v$class == "excluded", , drop = FALSE], n)
        if (nrow(ex)) {
            for (i in seq_len(nrow(ex))) cat(sprintf("    - %s [%s]\n", substr(ex$pathway[i], 1, 56), ex$reason[i]))
        } else cat("    (none)\n")
        cat("\n  undetermined:\n")
        u <- utils::head(v[v$class == "undetermined", , drop = FALSE], n)
        if (nrow(u)) {
            for (i in seq_len(nrow(u))) cat(sprintf("    - %s [%s]\n", substr(u$pathway[i], 1, 56), u$reason[i]))
        } else cat("    (none)\n")
    }
    if (!is.null(x$gaps) && nrow(x$gaps)) {
        cat("\n  driver inventory:\n")
        for (i in seq_len(nrow(x$gaps))) {
            g <- x$gaps[i, ]
            cat(sprintf("    %s: %d dataset(s), %d cell line(s)%s\n", g$driver,
                        g$n_datasets, g$n_cell_lines,
                        if (is.na(g$gap)) "" else paste0(" -- ", g$gap)))
        }
    }
    if (nrow(x$next_steps)) {
        cat("\n  next experiment:\n")
        for (i in seq_len(nrow(x$next_steps))) {
            cat(sprintf("    %d. %s [%s, unlocks %d]\n", i, x$next_steps$action[i],
                        x$next_steps$side[i], x$next_steps$unlocks[i]))
        }
    }
    invisible(x)
}
