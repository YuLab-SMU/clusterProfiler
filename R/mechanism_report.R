## Report evidence limits without turning missing support into causal exclusion.
.mechanism_remedies <- function() {
    rows <- list(
        list("single-cell-line", "data", "test the candidate in an independent, relevant cell context", TRUE),
        list("no-dataset", "data", "no usable dataset was resolved; check coverage or obtain a relevant perturbation", TRUE),
        list("lookup-failed", "data", "the lookup could not run (network or missing gpsadb); retry before reading these rows", TRUE),
        list("no-test-completed", "data", "the calls failed; rerun before interpreting the evidence", TRUE),
        list("too-small", "analysis", "check member coverage and whether a better-powered pathway test is appropriate", TRUE),
        list("low-support", "analysis", "check the contrast and pathway direction; do not lower thresholds merely to obtain support", FALSE),
        list("member-derived", "analysis", "state how pathway direction was defined; member agreement alone is not independent evidence", FALSE),
        list("no-members", "analysis", "check the identifier space and measured member coverage", FALSE),
        list("no-fold-changes", "analysis", "supply the measured fold changes to assess pathway movement", FALSE),
        list("no-direction", "analysis", "inspect mixed member directions before proposing a directional mechanism", FALSE),
        list("not-significant", "concept", "no significant support was detected in the tested datasets; assess power and context before drawing a negative conclusion", FALSE),
        list("coherent-refuted", "concept", "the legacy verdict means no significant support, not that the candidate was disproved", FALSE),
        list("no-perturbation-test", "analysis", "resolve a candidate and usable perturbation data before testing", FALSE),
        list("direction-conflict", "concept", "the tested perturbation conflicts with the proposed directional explanation; inspect context and alternative roles", FALSE),
        list("supported", "data", "validate the candidate in an independent perturbation; consistency is not causal identification", FALSE))
    data.frame(reason = vapply(rows, `[[`, character(1), 1L),
               side = vapply(rows, `[[`, character(1), 2L),
               action = vapply(rows, `[[`, character(1), 3L),
               unlocks = vapply(rows, `[[`, logical(1), 4L),
               stringsAsFactors = FALSE)
}

.mechanism_remedy_for <- function(reason, verdict) {
    key <- reason
    if (is.na(key) || !nzchar(key)) {
        key <- if (identical(verdict, "perturbation-supported")) "supported"
            else if (identical(verdict, "too-small")) "too-small" else NA_character_
    }
    tab <- .mechanism_remedies()
    hit <- tab[match(key, tab$reason), , drop = FALSE]
    if (!nrow(hit) || is.na(hit$side))
        return(list(side = NA_character_, action = NA_character_, unlocks = FALSE))
    list(side = hit$side, action = hit$action, unlocks = hit$unlocks)
}

.mechanism_verdict_class <- function(x) {
    if (!"level" %in% names(x)) return(rep("undetermined", nrow(x)))
    v <- if ("verdict" %in% names(x)) x$verdict else rep(NA_character_, nrow(x))
    ifelse(!is.na(x$level) & x$level == "A", "supported",
           ifelse(v %in% "direction-conflict", "conflicting", "undetermined"))
}

.mechanism_driver_names <- function(gene) {
    aliases <- tryCatch(mechanism_regulator_aliases(), error = function(e) character(0))
    back <- if (length(aliases)) names(aliases)[aliases == gene] else character(0)
    fwd <- tryCatch(mechanism_alias_regulator(gene), error = function(e) gene)
    tries <- unique(c(gene, fwd, back))
    tries[!is.na(tries) & nzchar(tries)]
}

.mechanism_inventory <- function(drivers, inventory = NULL) {
    query_status <- attr(inventory, "query_status")
    if (is.null(inventory)) {
        if (!length(drivers)) return(NULL)
        parts <- list()
        states <- vector("list", length(drivers))
        available <- requireNamespace("gpsadb", quietly = TRUE)
        for (i in seq_along(drivers)) {
            g <- drivers[i]
            found <- NULL
            errors <- character()
            if (available) for (nm in .mechanism_driver_names(g)) {
                d <- tryCatch(as.data.frame(
                    gpsadb::gpsadb_search_datasets(perturbed_gene = nm, limit = 100)[[1]]),
                    error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
                if (!is.null(d) && nrow(d) && "method" %in% names(d))
                    d <- d[!is.na(d$method) & !grepl("overexpress", d$method, ignore.case = TRUE), , drop = FALSE]
                if (!is.null(d) && nrow(d)) {
                    found <- d
                    found$queried_as <- g
                    break
                }
            }
            status <- if (!available) "unavailable" else if (!is.null(found)) "available"
                else if (length(errors)) "failed" else "empty"
            states[[i]] <- data.frame(driver = g, status = status,
                message = if (!available) "gpsadb is not installed" else paste(unique(errors), collapse = "; "),
                stringsAsFactors = FALSE)
            if (!is.null(found)) parts[[length(parts) + 1L]] <- found
        }
        query_status <- do.call(rbind, states)
        inventory <- if (length(parts)) do.call(rbind, parts) else
            data.frame(perturbed_gene = character(), dataset_id = character(),
                       cell_line = character(), stringsAsFactors = FALSE)
    }
    inventory <- as.data.frame(inventory, stringsAsFactors = FALSE)
    if (!"perturbed_gene" %in% names(inventory))
        rlang::abort("`inventory` needs a `perturbed_gene` column; a missing column is not evidence of absent datasets.")
    keep <- intersect(c("dataset_id", "perturbed_gene", "cell_line", "tissue", "method",
                        "n_samples", "in_reference", "queried_as", "n_up", "n_down"), names(inventory))
    inventory <- inventory[, keep, drop = FALSE]
    if ("dataset_id" %in% names(inventory)) {
        id <- as.character(inventory$dataset_id)
        known <- !is.na(id) & nzchar(trimws(id))
        inventory <- inventory[!known | !duplicated(id), , drop = FALSE]
    }
    attr(inventory, "query_status") <- query_status
    inventory
}

## These are unsigned, pairwise target-overlap neighbourhoods, not causal modules.
.mechanism_module_map <- function(drivers, network, min_overlap = 0.10) {
    if (!is.numeric(min_overlap) || length(min_overlap) != 1L || !is.finite(min_overlap) ||
        min_overlap < 0 || min_overlap > 1)
        stop("`min_overlap` must be one finite number between 0 and 1.", call. = FALSE)
    if (is.null(network) || !length(drivers)) return(NULL)
    net <- as.data.frame(network, stringsAsFactors = FALSE)
    if (!all(c("source", "target") %in% names(net)))
        rlang::abort("`network` needs `source` and `target` columns.")
    net$source <- as.character(net$source)
    net$target <- as.character(net$target)
    net <- net[!is.na(net$source) & nzchar(trimws(net$source)) &
               !is.na(net$target) & nzchar(trimws(net$target)), , drop = FALSE]
    targets <- lapply(split(net$target, net$source), unique)
    overlap <- function(a, b) {
        A <- targets[[a]]; B <- targets[[b]]
        if (is.null(A) || is.null(B)) return(NA_real_)
        u <- length(union(A, B))
        if (!u) return(NA_real_)
        length(intersect(A, B)) / u
    }
    out <- lapply(drivers, function(d) {
        others <- setdiff(drivers, d)
        j <- vapply(others, function(o) overlap(d, o), numeric(1))
        c(d, others[!is.na(j) & j >= min_overlap])
    })
    stats::setNames(out, drivers)
}

.mechanism_global_effect <- function(inventory, drivers, threshold = 3000L) {
    if (is.null(inventory) || !nrow(inventory) || !length(drivers) ||
        !all(c("n_up", "n_down") %in% names(inventory))) return(NULL)
    breadth <- vapply(drivers, function(d) {
        keep <- inventory$perturbed_gene %in% .mechanism_driver_names(d)
        if ("queried_as" %in% names(inventory)) keep <- keep | inventory$queried_as %in% d
        v <- inventory$n_up[keep] + inventory$n_down[keep]
        v <- v[is.finite(v) & v >= 0]
        if (length(v)) stats::median(v) else NA_real_
    }, numeric(1))
    stats::setNames(breadth, drivers)
}

.mechanism_global_effect_reason <- function(driver, breadth, threshold) {
    if (is.null(breadth) || !is.finite(breadth) || breadth <= threshold) return(NA_character_)
    sprintf("global-effect: available perturbations of %s changed a median of about %s genes; this is a breadth flag, not proof of essentiality or a nonspecific effect on this pathway. Counts depend on context and DE thresholds.",
            driver, format(as.integer(breadth), big.mark = ","))
}

.mechanism_gap_table <- function(inventory, drivers, context, modules = NULL) {
    if (!length(drivers)) return(NULL)
    tissue <- if (is.list(context)) context$tissue else NULL
    status <- attr(inventory, "query_status")
    rows <- lapply(drivers, function(g) {
        state <- if (!is.null(status)) status$status[match(g, status$driver)] else NA_character_
        unknown <- is.null(inventory) || (!is.na(state) && state %in% c("failed", "unavailable"))
        d <- if (is.null(inventory)) data.frame() else
            inventory[inventory$perturbed_gene %in% .mechanism_driver_names(g), , drop = FALSE]
        if (!nrow(d)) {
            gap <- if (unknown) "inventory coverage unknown; lookup unavailable or failed"
                else if (identical(state, "empty")) "no usable dataset resolved by this query"
                else "no matching dataset in the supplied inventory"
            return(data.frame(driver = g, n_datasets = if (unknown) NA_integer_ else 0L,
                n_cell_lines = if (unknown) NA_integer_ else 0L, tissues = NA_character_,
                in_context_tissue = NA_integer_, gap = gap, stringsAsFactors = FALSE))
        }
        cells <- if ("cell_line" %in% names(d)) trimws(as.character(d$cell_line)) else character()
        cells <- unique(cells[!is.na(cells) & nzchar(cells)])
        tissues <- if ("tissue" %in% names(d)) trimws(as.character(d$tissue)) else character()
        tissues <- unique(tissues[!is.na(tissues) & nzchar(tissues)])
        matched <- if (length(tissue) && length(tissues)) sum(tolower(tissues) == tolower(tissue)) else NA_integer_
        gap <- if (length(tissue) && !is.na(matched) && matched == 0L)
            sprintf("no available dataset from %s; available tissues: %s", tissue, paste(utils::head(tissues, 5), collapse = ", "))
            else if (!length(cells)) "cell-line identities unknown"
            else if (length(cells) < 2L) sprintf("one cell line only (%s)", cells) else NA_character_
        data.frame(driver = g, n_datasets = nrow(d), n_cell_lines = length(cells),
                   tissues = paste(utils::head(tissues, 5), collapse = ", "),
                   in_context_tissue = matched, gap = gap, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out$module <- unname(vapply(out$driver, function(g) {
        if (is.null(modules) || is.null(modules[[g]])) NA_character_ else paste(modules[[g]], collapse = ", ")
    }, character(1)))
    out$module_size <- unname(vapply(out$driver, function(g) {
        if (is.null(modules) || is.null(modules[[g]])) NA_integer_ else length(modules[[g]])
    }, integer(1)))
    out
}

.mechanism_candidate_evidence <- function(tests) {
    empty <- data.frame(pathway = character(), driver = character(), n_tests = integer(),
        n_support = integer(), n_conflict = integer(), n_failed = integer(),
        supporting_cells = character(), conflicting_cells = character(), assessment = character(),
        stringsAsFactors = FALSE)
    if (is.null(tests) || !nrow(tests)) return(empty)
    required <- c("pathway", "driver", "dataset", "cell_line", "status")
    if (!all(required %in% names(tests)))
        stop("`perturbation_tests` is missing required evidence columns.", call. = FALSE)
    groups <- unique(tests[, c("pathway", "driver"), drop = FALSE])
    rows <- lapply(seq_len(nrow(groups)), function(i) {
        z <- tests[tests$pathway %in% groups$pathway[i] & tests$driver %in% groups$driver[i], , drop = FALSE]
        cells <- function(status) {
            a <- trimws(as.character(z$cell_line[z$status %in% status]))
            paste(sort(unique(a[!is.na(a) & nzchar(a)])), collapse = ", ")
        }
        ns <- sum(z$status %in% "supports"); nc <- sum(z$status %in% "conflicts")
        nf <- sum(z$status %in% "failed")
        data.frame(pathway = groups$pathway[i], driver = groups$driver[i],
            n_tests = sum(!z$status %in% c("failed", "not-computable")),
            n_support = ns, n_conflict = nc, n_failed = nf,
            supporting_cells = cells("supports"), conflicting_cells = cells("conflicts"),
            assessment = if (ns && nc) "mixed" else if (ns) "support-only" else if (nc) "conflict-only"
                else if (any(z$status %in% "not-significant")) "no-significant-support" else "unavailable",
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' Report candidate consistency, conflicts, and evidence gaps
#'
#' Summarise measured evidence without treating a nonsignificant test as causal
#' exclusion. A supported row is consistency with a directional explanation,
#' not identification of a unique upstream driver. Significant conflicting tests
#' remain visible even when other datasets support that candidate.
#'
#' @param x A table from [pathway_evidence()]. Its `perturbation_tests` attribute
#'   retains all candidate/dataset outcomes, not just the best result per pathway.
#' @param drivers Drivers to inventory. Defaults to the drivers requested in `x`
#'   and those named in its rows or retained tests.
#' @param inventory Optional available-dataset table with `perturbed_gene` and
#'   optionally `dataset_id`, `cell_line`, `tissue`, `n_up` and `n_down`. Without
#'   it GPSAdb is queried (up to 100 rows per resolved spelling). A failed lookup
#'   is unknown coverage, not evidence that no datasets exist. A supplied table
#'   describes only its own scope, not an exhaustive database inventory.
#' @param context Named list. `context$tissue` is compared with available tissues.
#' @param preflight Optional [mechanism_preflight()] result. If the input screen
#'   does not support attribution, additional experiments are not automatically
#'   ranked; the report explains the current limit without claiming that future
#'   experiments could never change the evidence.
#' @param network Optional regulator-target table with `source` and `target`.
#'   `module` is an unsigned pairwise Jaccard neighbourhood among the requested
#'   candidates. It ignores edge weights and is not a disjoint clustering,
#'   proof of biological equivalence, or evidence that expression experiments
#'   cannot distinguish its members. No module accuracy gain is assumed.
#' @param min_overlap Jaccard threshold in [0, 1], default `0.10`. Lower thresholds
#'   enlarge neighbourhoods and their chance coverage; this is descriptive.
#' @param global_effect_de Median DE gene count above which perturbation breadth
#'   is flagged, default `3000`. Counts depend on tissue, method, sample size and
#'   DE thresholds. The flag does not prove essentiality or pathway nonspecificity,
#'   and does not change a grade or score.
#' @return A `mechanismReport` list with `verdicts` (classes `supported`,
#'   `conflicting`, `undetermined`), `gaps`, `next_steps`, `summary`, `context`,
#'   `inventory`, `perturbation_tests` and `candidate_evidence`. The last two
#'   expose all retained tests and counts of supporting/conflicting results per
#'   pathway/candidate. `next_steps$unlocks` is retained for compatibility and
#'   counts rows addressed by an action, not promised grade upgrades. The legacy
#'   `coherent-refuted` verdict is treated as undetermined, not disproved.
#' @seealso [pathway_evidence()], [mechanism_preflight()]
#' @export
mechanism_report <- function(x, drivers = NULL, inventory = NULL,
                             context = list(), preflight = NULL,
                             network = NULL, min_overlap = 0.10,
                             global_effect_de = 3000L) {
    if (!is.numeric(global_effect_de) || length(global_effect_de) != 1L ||
        !is.finite(global_effect_de) || global_effect_de < 0)
        stop("`global_effect_de` must be one nonnegative finite number.", call. = FALSE)
    if (!is.data.frame(x) || !nrow(x))
        return(.mechanism_report_empty("No pathway was scored; inspect the input and coverage before proposing a mechanism."))
    if (!all(c("level", "pathway") %in% names(x)))
        rlang::abort("`x` needs `level` and `pathway` columns from pathway_evidence().")
    tests <- attr(x, "perturbation_tests")
    if (!is.null(tests) && nrow(tests) && "pathway" %in% names(tests))
        tests <- tests[tests$pathway %in% x$pathway, , drop = FALSE]
    candidates <- .mechanism_candidate_evidence(tests)
    driver_of_row <- if ("driver" %in% names(x)) as.character(x$driver) else rep(NA_character_, nrow(x))
    if (is.null(drivers)) drivers <- unique(c(attr(x, "drivers"), driver_of_row, candidates$driver))
    drivers <- unique(as.character(drivers))
    drivers <- drivers[!is.na(drivers) & nzchar(trimws(drivers))]
    all_drivers <- unique(c(drivers, stats::na.omit(driver_of_row), candidates$driver))
    inv <- .mechanism_inventory(all_drivers, inventory)
    mod <- .mechanism_module_map(all_drivers, network, min_overlap)
    breadth <- .mechanism_global_effect(inv, all_drivers, global_effect_de)
    reason <- if ("reason" %in% names(x)) x$reason else rep(NA_character_, nrow(x))
    klass <- .mechanism_verdict_class(x)
    vcol <- if ("verdict" %in% names(x)) x$verdict else klass
    rem <- lapply(seq_len(nrow(x)), function(i) .mechanism_remedy_for(reason[i], vcol[i]))
    module_of <- function(d) {
        if (is.null(mod) || is.na(d) || is.null(mod[[d]])) NA_character_ else paste(mod[[d]], collapse = ", ")
    }
    module_size_of <- function(d) {
        if (is.null(mod) || is.na(d) || is.null(mod[[d]])) NA_integer_ else length(mod[[d]])
    }
    global_of <- function(d) {
        if (is.na(d) || is.null(breadth) || !d %in% names(breadth)) return(NA_character_)
        .mechanism_global_effect_reason(d, breadth[[d]], global_effect_de)
    }
    verdicts <- data.frame(pathway = x$pathway, class = klass, verdict = vcol,
        reason = reason, detail = driver_of_row,
        side = vapply(rem, `[[`, character(1), "side"),
        action = vapply(rem, `[[`, character(1), "action"),
        module = unname(vapply(driver_of_row, module_of, character(1))),
        module_size = unname(vapply(driver_of_row, module_size_of, integer(1))),
        global_effect_reason = unname(vapply(driver_of_row, global_of, character(1))),
        stringsAsFactors = FALSE)
    gaps <- .mechanism_gap_table(inv, drivers, context, mod)
    blocked <- !is.null(preflight) && inherits(preflight, "mechanismPreflight") &&
        !identical(preflight$tier, "attribution")
    actionable <- if (blocked) verdicts[0, , drop = FALSE] else
        verdicts[!is.na(verdicts$action) & !is.na(verdicts$side) & verdicts$side != "concept", , drop = FALSE]
    next_steps <- if (!nrow(actionable)) {
        data.frame(action = character(), side = character(), unlocks = integer(), examples = character(), stringsAsFactors = FALSE)
    } else {
        agg <- lapply(split(actionable, actionable$action), function(z)
            data.frame(action = z$action[1], side = z$side[1], unlocks = nrow(z),
                examples = paste(utils::head(z$pathway, 3), collapse = "; "), stringsAsFactors = FALSE))
        agg <- do.call(rbind, agg)
        agg[order(-agg$unlocks), , drop = FALSE]
    }
    n_mod <- sum(!is.na(verdicts$module_size) & verdicts$module_size > 1L)
    n_glob <- sum(!is.na(verdicts$global_effect_reason))
    summary <- c(
        if (blocked) sprintf("attribution is not supported by the current input screen (%s); this is a present evidence limit, not a claim that future experiments cannot help.", paste(sub(":.*", "", preflight$reasons), collapse = ", ")),
        sprintf("%d pathways: %d supported, %d conflicting, %d undetermined.", nrow(verdicts),
                sum(klass == "supported"), sum(klass == "conflicting"), sum(klass == "undetermined")),
        "Supported means consistency with tested perturbations, not identification of a unique causal driver.",
        if (n_mod) sprintf("%d row(s) have target-overlap neighbours; shared annotations do not establish indistinguishability or upstream order.", n_mod),
        if (n_glob) sprintf("%d row(s) have a perturbation breadth flag; inspect global_effect_reason and context before interpreting specificity.", n_glob),
        if (any(candidates$assessment == "mixed")) "Some candidate/pathway pairs have both supporting and conflicting tests; inspect candidate_evidence and perturbation_tests.",
        if (blocked) "Resolve the reported input or coverage limits before ranking additional perturbation experiments."
        else if (nrow(next_steps)) paste0("Next: ", next_steps$action[1], " (addresses ", next_steps$unlocks[1], " row(s); no upgrade is guaranteed).")
        else "No specific follow-up is ranked from the available evidence; this does not rule out useful future experiments.")
    structure(list(verdicts = verdicts, gaps = gaps, next_steps = next_steps,
        summary = summary, context = context, inventory = inv,
        perturbation_tests = tests, candidate_evidence = candidates), class = "mechanismReport")
}

.mechanism_report_empty <- function(reason) {
    structure(list(verdicts = data.frame(pathway = character(), class = character(),
        verdict = character(), reason = character(), detail = character(), side = character(),
        action = character(), module = character(), module_size = integer(), global_effect_reason = character()),
        gaps = NULL, next_steps = data.frame(action = character(), side = character(), unlocks = integer(), examples = character()),
        summary = reason, context = list(), inventory = NULL, perturbation_tests = NULL,
        candidate_evidence = .mechanism_candidate_evidence(NULL)), class = "mechanismReport")
}

#' @rdname mechanism_report
#' @param n Number of rows to show per section.
#' @param ... Ignored.
#' @export
print.mechanismReport <- function(x, n = 6L, ...) {
    cat("<mechanism report>\n")
    for (s in x$summary) cat("  ", s, "\n", sep = "")
    v <- x$verdicts
    if (nrow(v)) for (cl in c("supported", "conflicting", "undetermined")) {
        cat("\n  ", cl, ":\n", sep = "")
        z <- utils::head(v[v$class == cl, , drop = FALSE], n)
        if (!nrow(z)) cat("    (none)\n") else for (i in seq_len(nrow(z)))
            cat(sprintf("    - %s [%s] <- %s\n", substr(z$pathway[i], 1, 56), z$reason[i], z$detail[i]))
    }
    mod <- v[!is.na(v$module_size) & v$module_size > 1L, , drop = FALSE]
    if (nrow(mod)) {
        keys <- unique(vapply(strsplit(mod$module, ", ", fixed = TRUE), function(z) paste(sort(z), collapse = ", "), character(1)))
        cat("\n  target-overlap neighbourhoods (unsigned, not causal equivalence):\n")
        for (k in keys) cat(sprintf("    {%s}\n", k))
    }
    ge <- v[!is.na(v$global_effect_reason), , drop = FALSE]
    if (nrow(ge)) {
        ge <- ge[!duplicated(ge$detail), , drop = FALSE]
        cat("\n  perturbation breadth warnings:\n")
        for (i in seq_len(nrow(ge))) cat(sprintf("    %s: %s\n", ge$detail[i], sub("^global-effect: ", "", ge$global_effect_reason[i])))
    }
    cands <- x$candidate_evidence
    if (!is.null(cands) && nrow(cands)) {
        cat("\n  candidate evidence (support / conflict / failed):\n")
        for (i in seq_len(min(nrow(cands), n))) cat(sprintf("    %s <- %s: %d / %d / %d [%s]\n",
            cands$pathway[i], cands$driver[i], cands$n_support[i], cands$n_conflict[i], cands$n_failed[i], cands$assessment[i]))
        if (nrow(cands) > n) cat("    ... inspect $candidate_evidence and $perturbation_tests for all candidates\n")
    }
    if (!is.null(x$gaps) && nrow(x$gaps)) {
        cat("\n  driver inventory:\n")
        for (i in seq_len(nrow(x$gaps))) {
            g <- x$gaps[i, ]
            cat(sprintf("    %s: %s dataset(s), %s cell line(s)%s%s\n", g$driver,
                if (is.na(g$n_datasets)) "unknown" else as.character(g$n_datasets),
                if (is.na(g$n_cell_lines)) "unknown" else as.character(g$n_cell_lines),
                if (!is.na(g$module_size) && g$module_size > 1L) paste0(" [neighbours: ", g$module, "]") else "",
                if (is.na(g$gap)) "" else paste0(" -- ", g$gap)))
        }
    }
    if (nrow(x$next_steps)) {
        cat("\n  next experiment / analysis:\n")
        for (i in seq_len(nrow(x$next_steps))) cat(sprintf("    %d. %s [%s, addresses %d]\n",
            i, x$next_steps$action[i], x$next_steps$side[i], x$next_steps$unlocks[i]))
    }
    invisible(x)
}
