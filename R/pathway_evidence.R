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
    tok <- trimws(unlist(strsplit(as.character(result$geneID), "/")))
    ## 空 token 与 NA 都不是成员。不过滤的话：
    ##   * `"A/B/C/D/E/"` 会数出 6 个成员，`n_members` 与陈述里的数字都多 1
    ##     （实测带尾斜杠与空格时得到 6 而不是 5）；
    ##   * `geneID = NA` 会得到一个叫 "NA" 的成员，于是"没有成员"被报成 1 个。
    ## 两者都会影响 `n_members < 5` 这个 `too-small` 判据的边界。
    unique(tok[!is.na(tok) & nzchar(tok)])
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
    ## 先用包里已有的别名表归一，再查。
    ##
    ## 不归一会**静默查空**：GPSAdb 里 `GR` 有 0 个敲低数据集，`NR3C1` 有 1 个；
    ## 调用方写 `drivers = "GR"` 是自然的，得到的是"这个驱动者没有任何数据"，
    ## 输出与"数据里确实不支持它"长得一模一样。这一条在 airway 上的实测就是这样：
    ## 9 条通路全部 C，看不出是别名问题。
    ## 别名表是**单向**的（`NR3C1 -> GR`），而调用方与数据库可能各站一端：
    ## GPSAdb 里 `NR3C1` 有数据、`GR` 没有，而用户写 `drivers = "GR"` 最自然。
    ## 所以双向都试，谁先查到用谁。
    aliases <- tryCatch(mechanism_regulator_aliases(), error = function(e) character(0))
    back <- if (length(aliases)) names(aliases)[aliases == gene] else character(0)
    fwd <- tryCatch(mechanism_alias_regulator(gene), error = function(e) gene)
    tries <- unique(c(gene, fwd, back))
    tries <- tries[!is.na(tries) & nzchar(tries)]

    look <- function(g) {
        d <- tryCatch(as.data.frame(
            gpsadb::gpsadb_search_datasets(perturbed_gene = g, method = "knockdown",
                                           limit = 50)[[1]]),
            error = function(e) NULL)
        if (is.null(d) || !nrow(d)) NULL else d
    }
    d <- NULL
    for (g in tries) { d <- look(g); if (!is.null(d)) break }
    if (is.null(d)) return(NULL)
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
    if (!nrow(d)) {
        ## 空输入也要出声。没有这句时，用户拿到一个空 data.frame，无法分辨
        ## "富集本来就没结果"与"函数出问题了"——实测两个数据集（雌二醇、以及其他
        ## 富集不到任何通路的条件实验）就停在这一步。
        warning(paste("The enrichment result holds no pathways, so there is nothing to score.",
                      "That is usually the condition experiment failing to concentrate its response,",
                      "not a fault in pathway_evidence."), call. = FALSE)
        return(.pathway_evidence_empty())
    }
    ## 守卫放在 `head()` **之后**。只放在前面时，`max_pathways = 0L` 会得到 0 行的 d，
    ## 后面的 `do.call(rbind, list())` 返回 NULL，`attr(out, ...)` 就报
    ## "attempt to set an attribute on NULL"——一条指不出原因的错误。
    max_pathways <- max(1L, as.integer(max_pathways)[1])
    d <- utils::head(d, max_pathways)
    if (!nrow(d)) return(.pathway_evidence_empty())

    ## 标识空间：**自动识别并翻译**，不让调用方每次自己 bitr。
    ##
    ## ReactomePA 吐 Entrez、GPSAdb 只认 Symbol，而这两个来源的组合正是主路径。
    ## 早先要求调用方自己转，代价是有人忘了就得到一个"跑通了、什么都没发生"的结果
    ## ——0.8 秒、零次网络调用、全部 untestable。识别成本很低（看成员是否像纯数字），
    ## 让它自动做比让人记住更可靠。
    d <- .pathway_ids_to_symbols(d, signed)

    ## 标识空间守卫。
    ##
    ## ReactomePA 的 `geneID` 是 **Entrez**，而 GPSAdb 与多数差异表用的是 **Symbol**。
    ## 两者对不上时，每个通路的成员一个都匹配不上 `signed`，方向推不出来、支持率是 NA，
    ## 结果是**全部 untestable，但看起来像正常的空结论**——实测 0.8 秒跑完、一个网络
    ## 调用都没发。这种"跑通了、什么都没发生"比报错危险得多。
    ##
    ## 判据用"一条通路都没匹配上"：只要有一条匹配到，就说明空间是通的，
    ## 少数通路匹配不到是正常的（成员不在差异表里）。
    if (!is.null(signed) && length(signed) && nrow(d)) {
        any_hit <- any(vapply(seq_len(nrow(d)), function(i) {
            m <- .pathway_members(d[i, ])
            length(m) && any(m %in% names(signed))
        }, TRUE))
        if (!any_hit) {
            rlang::abort(c(
                "No pathway member matched the `signed` gene list.",
                x = "Every pathway came back untestable, which usually means the two are in different identifier spaces.",
                i = "ReactomePA returns Entrez IDs; pass symbols in both, or translate one side first."
            ))
        }
    }

    ## 候选数据集只解析一次，不按通路重复查——同一个驱动者在几十条通路上是同一批数据集。
    sets <- list()
    if (length(drivers)) {
        for (g in drivers) {
            ds <- .pathway_datasets(g, limit = max_datasets_per_driver)
            if (!is.null(ds)) sets[[g]] <- ds[, c("dataset_id", "cell_line"), drop = FALSE]
        }
    }

    ## A 级不可达这件事在**驱动者层面**就能确定，不必从每个通路的 `reason` 去反推。
    ##
    ## 早先的写法是等结果出来看原因分布，实测漏报：airway 那批通路里一半是方向冲突、
    ## 一半是成员分散，没有"单一主导原因"，于是明明 GR 只有一个细胞系却没报。
    ## 判据放在这里更直接、也不会误报。
    if (length(sets)) {
        cells_each <- vapply(sets, function(z) length(unique(stats::na.omit(z$cell_line))), 0L)
        if (all(cells_each < 2L)) {
            warning(sprintf(
                "Level A is unreachable: %s has knockdown data in a single cell line, and level A requires two. Results here cannot exceed B.",
                paste(sprintf("`%s` (%d)", names(cells_each), cells_each), collapse = ", ")),
                call. = FALSE)
        }
    } else if (length(drivers)) {
        warning(sprintf(
            "No knockdown data found in GPSAdb for %s; every pathway will be described without perturbation evidence.",
            paste0("`", paste(drivers, collapse = "`, `"), "`")), call. = FALSE)
    }

    ## 按 (数据集, 成员签名) 记忆化。
    ##
    ## 不能按通路缓存——测试结果取决于成员，换了通路就是另一次检验。但**不同通路
    ## 常常共享完全相同的成员集**：缺氧结果里 `response to hypoxia`、
    ## `response to decreased oxygen levels`、`response to oxygen levels` 三个
    ## 都是同一批 33 个基因，原来各测一遍。签名相同时复用是安全的，因为请求完全相同。
    memo <- new.env(parent = emptyenv())
    calls <- 0L
    fails <- 0L
    ask <- function(id, members) {
        key <- paste0(id, "|", paste(sort(members), collapse = ","))
        if (!is.null(memo[[key]])) return(memo[[key]])
        calls <<- calls + 1L
        r <- tryCatch(gpsadb::gpsadb_test_gene_set(id, members),
                      error = function(e) NULL)
        val <- if (is.null(r)) list(.failed = TRUE) else {
            t <- as.data.frame(r[[1]])
            if (!nrow(t)) list(.failed = TRUE) else
                list(padj = as.numeric(t$padj[1]), es = as.numeric(t$es[1]))
        }
        if (isTRUE(val$.failed)) fails <<- fails + 1L
        memo[[key]] <- val
        val
    }

    rows <- lapply(seq_len(nrow(d)), function(i) {
        members <- .pathway_members(d[i, ])
        dsrc <- .pathway_direction_source(d[i, ], signed)
        direction <- .pathway_direction(d[i, ], members, signed)
        support <- .pathway_support(members, signed, direction)
        desc <- d$Description[i] %||% d$ID[i] %||% NA_character_
        ## 先给默认值：下面的"无方向"分支会提前 return，那时检验还没跑，
        ## 但已经要填 `n_tests_run` 这一列。
        n_tested <- 0L

        ## 没有方向有两种，必须是两个结论，不是一个：
        ##   * `signed` 没给 —— 我们没有测过成员往哪边走，这是**没测**；
        ##   * `signed` 给了但成员对半分 —— 测过了，它们没有共同方向，这是**测了不成立**。
        ## 混成一句"没有方向"会把后者说成缺失，与 `absent` / `untested` 分开是同一条纪律。
        if (is.null(direction)) {
            untested <- is.null(signed) || !length(signed)
            ## 这里其实是**三个**结论，不是一个：
            ##   * 没给 `signed`            -> no-fold-changes（没测）
            ##   * 给了但一个成员都不在里面  -> no-members（这批基因不在差异表里）
            ##   * 给了、有成员、但正负对半  -> no-direction（测了，没有共同方向）
            ## 原来后两者共用 "Members split evenly"，于是一条 0 成员的通路
            ## （例如 `geneID` 为 NA 或全是不认识的名字）也说成"成员对半分"——
            ## 与 `why_c` 那边把 member-derived 并进 low-support 是同一类错误。
            hit <- if (untested) character(0) else members[members %in% names(signed)]
            why <- if (untested) "no-fold-changes"
                else if (!length(hit)) "no-members"
                else "no-direction"
            return(data.frame(
                pathway = desc, direction = NA_character_, n_members = length(members),
                support = support, driver = NA_character_, dataset = NA_character_,
                cell_line = NA_character_, es = NA_real_, padj = NA_real_,
                effect = NA_character_, n_tests_run = n_tested,
                n_cell_lines_supported = 0L, level = "C",
                verdict = if (untested) "untestable" else "unsupported",
                reason = why,
                statement = switch(why,
                    "no-fold-changes" = "No signed fold changes were supplied, so the pathway's movement is unknown.",
                    "no-members" = "None of this pathway's members appear in the signed gene list, so its movement is unknown.",
                    "no-direction" = "Members split evenly, so the pathway has no common direction to explain."),
                stringsAsFactors = FALSE))
        }

        ## 逐个候选驱动者的**每一份**数据集都检验，不只看最强的那一份。
        ##
        ## 原版只留 padj 最小的一个，于是 A 级等于"某一份数据支持"。实测这个含义太弱：
        ## 四个驱动者里三个换一个细胞系就不成立，而且效应量符号翻转
        ## （`results/stage3_refutation.tsv`）。所以要看得见全部数据集，才能判断复现。
        tested <- list()
        n_attempt <- 0L
        n_fail <- 0L
        for (g in names(sets)) {
            for (j in seq_len(nrow(sets[[g]]))) {
                id <- as.character(sets[[g]]$dataset_id[j])
                got <- ask(id, members)
                n_attempt <- n_attempt + 1L
                ## 调用失败与"调用成功但没结果"必须分开。
                ## 原来两者都是 `next`，于是一次断网会让 `tested` 变空，
                ## 再被解释成"成员方向不一致"——实测 HIF1A（4 个细胞系）在
                ## 模拟断网下得到 `low-support`、零报警，用户看不出是网络问题。
                if (isTRUE(got$.failed)) { n_fail <- n_fail + 1L; next }
                if (is.na(got$padj)) next
                ## 同一个数据集可能被两个驱动者指向同一个基因的别名，只记一次。
                already <- vapply(tested, function(z) identical(z$dataset, id) &&
                                  identical(z$gene, g), TRUE)
                if (any(already)) next
                tested[[length(tested) + 1L]] <- list(
                    gene = g, dataset = id,
                    cell_line = as.character(sets[[g]]$cell_line[j]),
                    es = got$es, padj = got$padj)
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
            level <- if (agrees && replicated) "A" else "B"
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
                effect = effect, n_tests_run = n_tested,
                n_cell_lines_supported = length(cells_agree),
                level = level, verdict = verdict,
                ## 降级原因要**写出来**。`single-dataset` 与 `direction-conflict` 都落在 B，
                ## 但一个是"证据不够多"，一个是"证据指向相反"——读者必须能分开。
                reason = if (level == "A") NA_character_
                    else if (!agrees) "direction-conflict" else "single-cell-line",
                statement = statement,
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
                effect = NA_character_, n_tests_run = n_tested,
                n_cell_lines_supported = 0L, level = "B", verdict = "coherent-untested",
                reason = "no-perturbation-test",
                statement = sprintf(
                    "%.0f%% of members moved %s, so the pathway moved; no perturbation test identified a driver.",
                    100 * support, direction),
                stringsAsFactors = FALSE))
        }

        ## 走到这里的 C 有三种完全不同的含义，必须分开写：
        ##   * 成员不在差异表里            -> no-members
        ##   * 成员对不上、方向无外部依据   -> low-support / member-derived
        ##   * 检验跑了但一份都没给出 p 值  -> too-small（通路太小）
        ##   * 检验跑了、都不显著          -> not-significant
        ## 原来它们全是 `C unsupported`，读者无法分辨"没找到"与"没法找"。
        ## 顺序有讲究：成员只有个位数时，**这才是要报的原因**。
        ## 一个 2 成员的通路，方向只能从这 2 个基因自己推出来，于是必然落到
        ## `low-support`——但读者看到"成员方向不一致"会去改方向策略，
        ## 实际问题是他给的差异基因太少。实测 MCF-7（23 个上调）与 A549 依托泊苷
        ## （16 个上调）都因此被误标成 low-support。
        ## `member-derived` 与 `low-support` **必须分开**。原来合成一个 "low-support"，
        ## 于是 support = 100% 的通路也印出"低于 50% 的下限"——上一轮我以为修了，
        ## 实际加的那条 `else` 分支永远走不到，因为 why_c 从不会取到它。
        ## 这次每个值都有明确来源，`statement` 用 `switch` 一一对应，不留兜底。
        ## `no-test-completed` 排在 `member-derived` 之前：一次调用全崩时，
        ## "检验没能跑起来"比"方向由成员推导"更该被看见——否则用户会以为自己
        ## 的输入有问题，实际是网络。实测断网时逐条 reason 曾显示 `member-derived`，
        ## 只有全局 warning 说得出真相，两层信息不一致。
        why_c <- if (is.na(support)) "no-members"
            else if (n_attempt > 0L && n_fail == n_attempt) "no-test-completed"
            else if (length(members) < 5L) "too-small"
            else if (identical(dsrc, "members")) "member-derived"
            else if (support < min_support) "low-support"
            else if (!length(tested)) {
                if (length(sets)) "too-small" else "no-dataset"
            } else "not-significant"

        data.frame(
            pathway = desc, direction = direction, n_members = length(members),
            support = support, driver = NA_character_, dataset = NA_character_,
            cell_line = NA_character_, es = NA_real_, padj = NA_real_,
            effect = NA_character_, n_tests_run = n_tested,
            n_cell_lines_supported = 0L, level = "C", verdict = "unsupported",
            reason = why_c,
            ## 陈述必须与 `why_c` 一致。原先只有一句"低于支持下限"，但 `support` 是 100%
            ## 时也走这里（原因是方向由成员自身推导），于是输出"100% 的成员上调，
            ## 低于 50% 的下限"这种自相矛盾的话。实测 airway 上就这么印出来了。
            statement = switch(why_c,
                "no-members" = "Members were not found in the signed gene list; the pathway's movement is unknown.",
                "too-small" = sprintf("Only %d members reached this pathway, too few to establish its direction or test a driver; supply more differentially expressed genes.",
                                      length(members)),
                "member-derived" = "The pathway's direction was inferred from these same members, so majority agreement cannot serve as evidence for it.",
                "low-support" = sprintf("%.0f%% of members moved %s, below the support floor of %.0f%%.",
                                        100 * support, direction, 100 * min_support),
                "no-test-completed" = "Every perturbation call for this pathway failed, so the absence of evidence may be an outage rather than a result.",
                "no-dataset" = "No knockdown dataset was available for the candidate drivers, so the pathway was described but not tested.",
                "not-significant" = sprintf("The drivers were tested against these %d members and none reached significance, so no driver is supported here.",
                                            length(members)),
                "unsupported"),
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, rows)

    ## 跨源一致这一终点已移除，原因写在设计文档 §6.4：
    ## 三种第二来源（PPI 密度、成员重叠、GO↔Reactome 官方映射）实测都不成立，
    ## 而这一列在主路径下恒为 `no-crossref`——它不影响任何判定，留着会让读者
    ## 以为做过一次验证。**一个恒定的空结论比没有这一列更误导。**
    ## 要恢复，先有能成立的第二来源，而不是先有这一列。

    attr(out, "network_calls") <- calls
    out <- out[order(match(out$level, c("A", "B", "C")), -out$support), , drop = FALSE]

    ## 结构性上限要**说出来**，不能只体现在等级上。
    ##
    ## 没有这一段时，`全部 B` / `全部 C` 与"测了、不支持"完全同形——而实测九次里
    ## 有五次的低等级来自外部缺数据（GR 只有 1 个细胞系、通路太小、差异基因太少），
    ## 不是方法判断的结果。读者看不出差别，就会把"没法看"读成"没有证据"。
    ## 全部调用失败时直接说，不看 reason 分布——这是一次调用层面的故障，
    ## 任何等级都不该被当真。
    outage <- calls > 0L && fails == calls
    if (outage) {
        warning(sprintf("All %d GPSAdb calls failed. No perturbation evidence was collected; every grade below reflects an outage, not the data. Check the connection and rerun.", calls),
                call. = FALSE)
    }
    ## 断网时不再重复第二条：`_pathway_warn_capped` 那句说的正是同一件事。
    ## 一次故障发两条同义 warning 只会让人以为是两个问题。
    if (!outage) .pathway_warn_capped(out, drivers)
    ## 必须设 class，否则 `print.pathway_evidence` 永远不会被派发。
    ## 实测这份代码从写出来到现在**一次都没执行过**：返回的是裸 data.frame，
    ## 用户看到的是 20 列的原始表格，`reason` / `capped by` 那些行一个都没印。
    ## 又一个"看起来写好了、其实没跑"的例子。
    class(out) <- c("pathway_evidence", "data.frame")
    out
}

.pathway_evidence_empty <- function() {
    out <- data.frame(
        pathway = character(), direction = character(), n_members = integer(),
        support = numeric(), driver = character(), dataset = character(),
        cell_line = character(), es = numeric(), padj = numeric(),
        effect = character(), n_tests_run = integer(),
        n_cell_lines_supported = integer(),
        level = character(), verdict = character(), reason = character(),
        statement = character(), stringsAsFactors = FALSE)
    ## 与正常返回保持同一 class。不设的话空结果不走 print.pathway_evidence，
    ## 用户在同一函数上看到两种输出形态。
    class(out) <- c("pathway_evidence", "data.frame")
    out
}

## 列子集后必须退回普通 data.frame。
##
## `[.data.frame` 会把自定义 class 一起带过去，于是 `out[, c("pathway","level")]`
## 仍然声称自己是 `pathway_evidence`，而 `print.pathway_evidence` 要读 `driver`、
## `reason`、`support` 这些已经不在的列，`is.na(NULL)` 是 `logical(0)`，
## `if` 收到长度 0 就报 "argument is of length zero"——一条指不出原因的错误。
##
## 按列取子集是极常见的操作（本项目自己的脚本里就到处都是），所以把 class 的
## 存续条件绑在"列还全"上，而不是让 print 去猜。
#' @export
`[.pathway_evidence` <- function(x, i, j, drop = TRUE) {
    needed <- c("pathway", "direction", "n_members", "support", "driver",
                "cell_line", "es", "padj", "level", "reason", "statement")
    out <- NextMethod()
    if (is.data.frame(out) && !all(needed %in% names(out))) {
        class(out) <- "data.frame"
    }
    out
}

#' @rdname pathway_evidence
#' @export
print.pathway_evidence <- function(x, n = 10L, ...) {
    cat("<pathway evidence>\n")
    counts <- table(factor(x$level, levels = c("A", "B", "C")))
    cat("  levels:       ", paste0(names(counts), " ", as.integer(counts),
                                  collapse = "  "), "\n", sep = "")

    ## 降级原因分行列出。等级只是结论，`reason` 才是"为什么停在这里"。
    rs <- x$reason[!is.na(x$reason)]
    if (length(rs)) {
        rc <- sort(table(rs), decreasing = TRUE)
        cat("  capped by:    ", paste0(names(rc), " ", as.integer(rc),
                                       collapse = ", "), "\n", sep = "")
    }

    shown <- utils::head(x, n)
    if (nrow(shown)) {
        cat("\n")
        for (i in seq_len(nrow(shown))) {
            row <- shown[i, ]
            driver <- if (!length(row$driver) || is.na(row$driver)) "-" else
                sprintf("%s/%s padj %.2g", row$driver, row$cell_line, row$padj)
            ## `%||%` 只接 NULL，接不住 NA；原来方向印成 "NA"、支撑率印成 "NA%"。
            dir_txt <- if (is.na(row$direction)) "?" else row$direction
            sup_txt <- if (is.na(row$support)) "-" else sprintf("%.0f%%", 100 * row$support)
            cat(sprintf("  [%s] %s (%s, %d members, %s together)\n",
                        row$level, substr(row$pathway, 1, 52), dir_txt,
                        row$n_members, sup_txt))
            why <- if (length(row$reason) && !is.na(row$reason)) sprintf("  [%s]", row$reason) else ""
            cat(sprintf("      driver: %s%s\n", driver, why))
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


## 把通路结果里的成员 ID 统一到符号空间。
##
## 判定：取前若干通路的成员，若**全部**是纯数字，就认定是 Entrez 并翻译。
## 只做这一种判定，因为 Entrez/RefSeq/Ensembl 之外没有常见的纯数字形态，
## 而猜错的代价（静默全空）比不猜的代价（要求调用方自己转）大得多。
##
## 翻译失败时**原样返回**，让后面的守卫去报错——错误信息在那里更完整。
.pathway_ids_to_symbols <- function(d, signed) {
    if (is.null(d) || !nrow(d) || !"geneID" %in% names(d)) return(d)
    probe <- utils::head(d$geneID, 5L)
    toks <- unique(unlist(strsplit(as.character(probe), "/")))
    toks <- toks[nzchar(toks) & !is.na(toks)]
    if (!length(toks) || !all(grepl("^[0-9]+$", toks))) return(d)
    ## 已经是数字、且 signed 也是数字时不用转，两者本来就同空间。
    if (!is.null(signed) && length(signed) && all(grepl("^[0-9]+$", names(signed)))) return(d)
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(d)

    conv <- tryCatch(suppressMessages(
        AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = toks,
                              keytype = "ENTREZID", columns = "SYMBOL")),
        error = function(e) NULL)
    if (is.null(conv) || !nrow(conv)) return(d)
    map <- stats::setNames(as.character(conv$SYMBOL), as.character(conv$ENTREZID))
    d$geneID <- vapply(strsplit(as.character(d$geneID), "/"), function(ids) {
        sym <- unname(map[ids])
        paste(stats::na.omit(sym), collapse = "/")
    }, character(1))
    d[nzchar(d$geneID), , drop = FALSE]
}


## 当一组通路的低等级**全部**来自同一个外部缺口时，发一次 warning 说清楚。
##
## 只在"没有任何 A"时触发：如果已经有 A，说明条件部分满足，逐行看 `reason` 就够了，
## 再报一次只是噪音。
.pathway_warn_capped <- function(out, drivers) {
    ## `!any(level %in% c("A","B","C"))` 是恒假的：level 只有这三个取值。
    ## 原来那句永远不会命中，等于没写。
    if (!nrow(out)) return(invisible(NULL))
    if (any(out$level == "A", na.rm = TRUE)) return(invisible(NULL))

    why <- c("single-cell-line", "no-dataset", "too-small", "no-test-completed")
    rs <- out$reason[!is.na(out$reason) & out$reason %in% why]
    if (!length(rs) || length(unique(rs)) != 1L) return(invisible(NULL))
    ## 要占多数才报——零星几条小通路不够格解释整批结果。
    if (length(rs) < max(1L, floor(0.5 * nrow(out)))) return(invisible(NULL))

    who <- paste0("`", paste(drivers, collapse = "`, `"), "`")
    msg <- switch(unique(rs),
        "single-cell-line" = sprintf(
            "level A is unreachable: %s has knockdown data in one cell line, and level A needs two. Everything is capped at B.",
            who),
        "no-dataset" = sprintf(
            "level A is unreachable: GPSAdb returned no knockdown dataset for %s. The pathways are described, not tested.",
            who),
        "too-small" = paste(
            "level A is unreachable: the pathways were too small for test_gene_set to produce a p-value,",
            "so the perturbation evidence could not be evaluated at all."),
        "no-test-completed" = paste(
            "every GPSAdb call failed, so no perturbation evidence could be collected at all.",
            "The C grades here describe an outage, not the biology -- check the connection and rerun.")
    )
    warning(msg, call. = FALSE)
    invisible(NULL)
}
