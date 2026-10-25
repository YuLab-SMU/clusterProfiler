#' Collapse redundant findings before reporting them
#'
#' Groups findings whose gene sets largely overlap and keeps the strongest
#' member of each group, so a result reads as a handful of statements rather
#' than as every set that crossed a threshold.
#'
#' The need is not uniform across sources, which is why this is a separate step
#' rather than something [consensus()] always does. Measured on the airway list:
#' transcription-factor findings barely overlap at all (median pairwise Jaccard
#' 0.000 over 155 findings), cell-state findings likewise (0.008 over 337), and
#' only tissue signatures are genuinely redundant (17% of pairs above 0.5, and
#' collapsing them keeps 17 of 42). Running a collapse step that is not needed
#' discards findings that say different things.
#'
#' @param x An `evidenceBundle` from [evidence()].
#' @param threshold Jaccard index above which two findings are the same
#'   statement. 0.5 means "shares most of its genes".
#' @param species Species whose collections supply the gene sets.
#' @param release MSigDB release, taken from the bundle when `NULL`.
#' @return A data frame like `consensus_findings()`, with two extra columns:
#'   `collapsed` (how many findings this row stands for) and `absorbed` (the
#'   labels it replaced). Ordered by significance.
#' @seealso [consensus()], [consensus_findings()]
#' @examples
#' \dontrun{
#' ev <- evidence(genes, upstream = "tissue")
#' collapse_findings(ev)
#' }
#' @export
collapse_findings <- function(x,
                              threshold = 0.5,
                              species = NULL,
                              release = NULL) {
    if (!inherits(x, "evidenceBundle")) {
        rlang::abort("`x` must be an evidenceBundle from evidence().")
    }
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold <= 0 || threshold > 1) {
        rlang::abort("`threshold` must be one number in (0, 1].")
    }

    findings <- consensus_findings(consensus(x))
    if (!nrow(findings)) {
        return(.collapse_empty())
    }

    if (is.null(species)) {
        species <- x$input$design$species %||% "Homo sapiens"
    }
    if (is.null(release)) {
        release <- unique(stats::na.omit(x$provenance$release))
        release <- if (length(release)) release[1] else "2026.2.Hs"
    }

    ## Gene sets come back keyed by collection term; a finding is the union of
    ## its terms, so the comparison happens at the level the verdict was made.
    by_term <- .collapse_members(x, findings, species, release)
    if (!length(by_term)) return(.collapse_empty())
    term_of <- .collapse_finding_terms(x, findings)
    members <- lapply(term_of, function(terms) {
        genes <- unique(unlist(by_term[intersect(terms, names(by_term))]))
        if (length(genes)) genes else NULL
    })
    members <- members[!vapply(members, is.null, TRUE)]

    keep <- .collapse_greedy(findings, members, threshold)

    ## `keep` holds row indices into `findings`; the absorbed findings were
    ## attributed during the greedy pass, so nothing is counted twice.
    out <- findings[keep, , drop = FALSE]
    taken <- attr(keep, "absorbed") %||% list()
    out$collapsed <- vapply(seq_along(keep), function(j) {
        1L + length(taken[[as.character(keep[j])]] %||% character(0))
    }, integer(1))
    out$absorbed <- vapply(seq_along(keep), function(j) {
        got <- taken[[as.character(keep[j])]] %||% character(0)
        if (!length(got)) NA_character_ else paste(got, collapse = ", ")
    }, character(1))
    rownames(out) <- NULL
    out
}

#' @rdname collapse_findings
#' @param x An object to print.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.collapsedFindings <- function(x, ...) {
    cat("<collapsedFindings>\n")
    cat("  kept:     ", nrow(x), " finding(s)\n", sep = "")
    n_absorbed <- sum(x$collapsed > 1L)
    if (n_absorbed) {
        cat("  absorbed: ", sum(x$collapsed) - nrow(x),
            " redundant finding(s)\n", sep = "")
    }
    if (nrow(x)) {
        cat("\n")
        show <- utils::head(x, 10)
        print(
            data.frame(
                regulator = show$regulator,
                verdict   = show$verdict,
                level     = show$max_level,
                padj      = signif(show$padj, 3),
                collapsed = show$collapsed,
                stringsAsFactors = FALSE
            ),
            row.names = FALSE
        )
    }
    invisible(x)
}

## The empty result, with exactly the columns a populated one has.
##
## Written out rather than derived, so this is the one place that can drift
## from consensus(). `evidence_kind` was added to the verdict table and missed
## here, which made an empty result and a non-empty one disagree about their
## own shape -- the kind of difference that only shows up in a caller that
## happens to get no findings.
.collapse_empty <- function() {
    out <- data.frame(
        regulator = character(), sources = character(), n_sources = integer(),
        max_level = character(), verdict = character(), padj = numeric(),
        family = logical(), evidence_kind = character(),
        direction = character(), writable = character(),
        collapsed = integer(), absorbed = character(), stringsAsFactors = FALSE)
    structure(out, class = c("collapsedFindings", "data.frame"))
}

## Gene membership for each finding, looked up from the collections.
##
## The bundle stores only collection metadata, not the gene lists, so the sets
## are rebuilt. They come from the same cache evidence() used, so a cached run
## pays no download.
.collapse_members <- function(x, findings, species, release) {
    ## The gene set is keyed by the collection term, not by the regulator label:
    ## a finding's `regulator` is the collapsed entity name ("AP1") while the
    ## set it came from is the original term ("TGANTCA_AP1_C"). Indexing on the
    ## label matched only the minority of findings whose label happens to equal
    ## its term, and silently dropped the rest.
    r <- x$regulators
    r <- r[r$regulator %in% findings$regulator, , drop = FALSE]
    r <- r[!is.na(r$padj_within_source) & r$padj_within_source < 0.05, , drop = FALSE]
    sources <- intersect(unique(r$source), c("tf", "mirna", "cellstate", "tissue"))
    if (!length(sources) || !nrow(r)) return(list())

    out <- list()
    for (src in sources) {
        g <- tryCatch(mechanism_sets(src, species = species, release = release),
                      error = function(e) NULL)
        if (is.null(g)) next
        m <- g@gsid2gene
        for (term in unique(r$term[r$source == src])) {
            genes <- unique(m$gene[m$gsid == term])
            if (length(genes)) out[[term]] <- genes
        }
    }
    out
}

## Map each finding to the terms that produced it.
##
## A finding can come from several terms (AP1 is enriched through three motif
## sets), and collapsing compares the finding as a whole, so its membership is
## the union over those terms.
.collapse_finding_terms <- function(x, findings) {
    r <- x$regulators
    r <- r[r$regulator %in% findings$regulator, , drop = FALSE]
    r <- r[!is.na(r$padj_within_source) & r$padj_within_source < 0.05, , drop = FALSE]
    split(r$term, r$regulator)
}

## Greedy collapse in significance order: keep a finding unless it restates one
## already kept. Greedy rather than clustering because the output must be
## ordered by strength, and the strongest member of a group is the one worth
## reporting.
.collapse_greedy <- function(findings, members, threshold) {
    ## `members` is keyed by finding label, built from the union of that
    ## finding's terms; a finding with no gene set at all is still kept, since
    ## there is nothing to show it restates.
    ordered <- order(findings$padj, na.last = TRUE)
    if (!length(ordered)) return(integer(0))

    keep <- integer(0)
    ## Each dropped finding is attributed to the kept finding that absorbed it,
    ## decided at the moment it is dropped. Recomputing "how many were similar"
    ## afterwards double-counts, because two kept findings can both be similar
    ## to the same dropped one.
    absorbed_by <- list()
    for (i in ordered) {
        mine <- members[[findings$regulator[i]]]
        if (is.null(mine)) {
            keep <- c(keep, i)
            next
        }
        hit <- NA_integer_
        for (k in keep) {
            other <- members[[findings$regulator[k]]]
            if (!is.null(other) && .collapse_jaccard(mine, other) > threshold) {
                hit <- k
                break
            }
        }
        if (is.na(hit)) {
            keep <- c(keep, i)
        } else {
            key <- as.character(hit)
            absorbed_by[[key]] <- c(absorbed_by[[key]],
                                    findings$regulator[i])
        }
    }
    attr(keep, "absorbed") <- absorbed_by
    keep
}

.collapse_jaccard <- function(a, b) {
    if (!length(a) || !length(b)) return(0)
    inter <- length(intersect(a, b))
    if (!inter) return(0)
    inter / length(union(a, b))
}
