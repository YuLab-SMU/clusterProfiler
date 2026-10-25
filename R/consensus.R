#' Decide which mechanisms the evidence actually supports
#'
#' Combine an `evidenceBundle` into one verdict per regulator. The point is not
#' to rank by p-value but to say which claims survive a change of evidence
#' type, and -- just as importantly -- which ones were never examined.
#'
#' Four verdicts:
#'
#' * `supported` -- at least `min_independent` distinct sources agree **on the
#'   same entity**, and at least one of them is graded B or above. A regulator
#'   found only at level C is never `supported`, however small its p-value.
#'
#'   Agreement is only counted for an identical regulator label. A transcription
#'   factor and a microRNA that happen to share targets are two different
#'   entities: `FOXO3` supported by TF target sets and `MIR23C` supported by
#'   miRNA target sets are each `single_source`, not one `supported` finding.
#'   Cross-type agreement is a statement about the gene list, not about either
#'   regulator, and it belongs in [evidence_ppi()] or a manual comparison.
#'
#'   The same applies across levels of the same source family: a cell-state name
#'   such as `HAY_BONE_MARROW_STROMAL` and a tissue name such as `CAR_HPX` are
#'   entities to be compared with **each other**, never summed with a factor.
#'   Their identifiers also mean something different -- "this signature
#'   resembles that state" rather than "this factor may drive these genes" --
#'   so a `supported` verdict for them would be a category error, not evidence.
#'   Both are graded C by [mechanism_xref()] for that reason, and a C-level
#'   entity is reported as single_source at best.
#'
#'   **This has a consequence worth stating plainly.** With `upstream = c("tf",
#'   "mirna")` two sources never name the same entity, so `supported` cannot
#'   occur; every finding is `single_source`, `absent` or `untested`. The verdict
#'   becomes reachable once evidence that names the *same* entity is added --
#'   a second target collection for the same factor, a perturbation signature
#'   naming the feature it perturbed, or a network grade attached to the same
#'   regulator. Until then `consensus()` still does its job: it separates what
#'   was tested and not found from what was never tested, which is the part a
#'   p-value threshold cannot do.
#' * `single_source` -- one source names it. That is a hypothesis, not a
#'   mechanism.
#' * `absent` -- the source could have spoken about these genes (its coverage
#'   reaches `coverage_floor`, or no floor was set) and did not.
#' * `untested` -- the source never examined these genes. **`absent` is
#'   evidence; `untested` is the absence of evidence**, and reporting the second
#'   as the first is the commonest way a regulator analysis overstates itself.
#'
#' @param x An `evidenceBundle` from [evidence()].
#' @param min_independent Distinct sources required for `supported`.
#' @param coverage_floor A source whose `coverage_tested` is below this is
#'   reported as `untested` whatever it found. `NULL` takes the floor recorded
#'   in `x$input$design`.
#' @return A `data.frame` with `regulator`, `sources`, `n_sources`,
#'   `max_level`, `verdict`, `padj`, `direction` and `writable` (the sentence
#'   the evidence permits; see [mechanism_writable()]).
#' @seealso [evidence()], [mechanism_writable()]
#' @examples
#' \dontrun{
#' ev <- evidence(c("HIF1A", "PGK1", "CA9"))
#' consensus(ev)
#' }
#' @export
consensus <- function(x, min_independent = 2, coverage_floor = NULL) {
    if (!inherits(x, "evidenceBundle")) {
        rlang::abort("`x` must be an evidenceBundle from evidence().")
    }
    if (is.null(coverage_floor)) coverage_floor <- x$input$design$coverage_floor

    r <- x$regulators
    if (!nrow(r)) {
        return(data.frame(
            regulator = character(), sources = character(), n_sources = integer(),
            max_level = character(), verdict = character(), padj = numeric(),
            family = logical(), evidence_kind = character(),
            direction = character(), writable = character(),
            stringsAsFactors = FALSE
        ))
    }

    ## A source that never examined these genes cannot report an absence.
    counts <- x$diagnostics$counts
    untested_source <- counts$source[
        !is.na(counts$coverage_floor) & counts$coverage_tested < counts$coverage_floor
    ]

    significant <- r[!is.na(r$padj_within_source) & r$padj_within_source < 0.05, , drop = FALSE]
    significant <- significant[!significant$source %in% untested_source, , drop = FALSE]
    ## A set that names no factor cannot support a claim about one: a consensus
    ## sequence is not a regulator, however strongly it is enriched.
    significant <- significant[significant$named, , drop = FALSE]

    regulators <- sort(unique(r$regulator))
    rows <- lapply(regulators, function(one) {
        mine <- significant[significant$regulator == one, , drop = FALSE]
        sources <- sort(unique(mine$source))
        n <- length(sources)

        ## "A" < "B" < "C" alphabetically, so min() is the strongest grade.
        ## na.rm keeps one ungraded set from erasing the whole regulator's
        ## level; if every level is missing the answer is NA, which
        ## mechanism_writable() renders as no sentence.
        levels <- mine$level[!is.na(mine$level)]
        max_level <- if (length(levels)) min(levels) else NA_character_

        ## A regulator that was tested and did not reach significance is
        ## absent, not untested: the source looked at it and said no. Untested
        ## is reserved for a regulator no live source tested at all -- which
        ## includes the case where the only source naming it was excluded for
        ## insufficient coverage.
        examined <- r[r$regulator == one & r$named & !r$source %in% untested_source, , drop = FALSE]

        ## Two sources agreeing on the string is not two sources agreeing on
        ## the entity: TRANSFAC's NF1 is Nuclear Factor I and the HGNC symbol
        ## NF1 is neurofibromin, so a motif set and a knockout dataset "agree"
        ## while describing different molecules. A collided label therefore
        ## cannot reach `supported`.
        collided <- mechanism_is_collision(one)

        ## What makes corroboration informative is not how many sources agree
        ## but what they measured. Three rules, each from a failure mode:
        ##
        ## * At least one row must be `measured`. Two predictions agreeing
        ##   (motif + miRDB) is one observation confirmed by the same reasoning
        ##   twice, and a prediction plus a tissue signature is two things that
        ##   were never perturbed. Neither is evidence that a factor drives the
        ##   list.
        ## * Independent measurements corroborate even when they share a kind.
        ##   Two knockdown datasets are the strongest position available, so
        ##   counting distinct kinds instead of distinct sources under-reports
        ##   exactly the best case.
        ## * Observation alone never corroborates: a signature says the genes
        ##   resemble a tissue, not that anything regulates them.
        kinds <- mine$evidence_kind %||% character(0)
        kinds <- kinds[!is.na(kinds) & nzchar(kinds)]
        ## "none" is an unnamed motif: it names no factor, so it cannot count.
        kinds <- setdiff(kinds, "none")
        n_measured <- sum(kinds == "measured")

        ## A source counts once. Two motif sets in the same collection are one
        ## witness, and two rows from one source are one source.
        n_witness <- n

        corroborated <- n_witness >= min_independent &&
            n_measured >= 1L &&
            n_measured >= min(n_witness, 1L) &&
            any(mine$level %in% c("A", "B"))

        verdict <- if (!nrow(examined)) {
            "untested"
        } else if (!n) {
            "absent"
        } else if (!collided && corroborated) {
            "supported"
        } else if (n >= 1) {
            "single_source"
        } else {
            "untested"
        }

        data.frame(
            regulator  = one,
            sources    = paste(sources, collapse = ","),
            n_sources  = n,
            max_level  = max_level,
            verdict    = verdict,
            padj       = if (n) min(mine$padj_within_source) else NA_real_,
            ## Inherited from the evidence rows: TRUE when this label names a
            ## family or complex, which no perturbation dataset can cover.
            family     = isTRUE(any(examined$family %||% FALSE)),
            ## The kinds that produced this verdict, so a reader can see whether
            ## corroboration came from different kinds of evidence or from the
            ## same kind twice.
            ## Distinct kinds, so "measured,measured" reads as one word: the
            ## count of witnesses is already `n_sources`.
            evidence_kind = if (length(kinds)) {
                paste(sort(unique(kinds)), collapse = ",")
            } else NA_character_,
            ## Only a measured perturbation states a direction. Where several
            ## datasets disagree the entity is "mixed", and a factor with no
            ## measurement keeps NA rather than guessing from sign of a score.
            direction  = {
                dirs <- unique(mine$direction[!is.na(mine$direction)])
                if (!length(dirs)) NA_character_
                else if (length(dirs) == 1L) dirs
                else "mixed"
            },
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out <- out[order(out$padj, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL
    out$writable <- mechanism_writable(out$verdict, out$max_level, out$regulator,
                                       sources = out$sources, family = out$family)
    out
}

## The subset of a verdict table worth reading: entities the evidence actually
## speaks about. An `absent` verdict is reported once per source in
## `$diagnostics`; repeating it for every one of a thousand factor labels turns
## the sentence table into noise and hides the few that matter.
consensus_findings <- function(cv) {
    keep <- cv$verdict %in% c("supported", "single_source")
    cv[keep, , drop = FALSE]
}

#' The sentence an evidence level permits
#'
#' Turn a verdict and a level into the strongest statement the evidence
#' supports, with the regulator name filled in. This is what keeps a level C
#' finding out of a mechanism sentence: the wording available at C describes a
#' correlation, and no higher wording is offered.
#'
#' @param verdict A `consensus()` verdict.
#' @param level An evidence level (`"A"`, `"B"`, `"C"`) or `NA`.
#' @param regulator The entity name.
#' @param sources Comma-separated source labels as `consensus()` records them.
#'   A measured perturbation is graded A whether it hit a regulator or a
#'   downstream effector, so the sentence for a single-source A finding depends
#'   on whether a target-set collection also names the gene as a regulator.
#' @return A character vector of sentences.
#' @seealso [consensus()]
#' @examples
#' mechanism_writable(c("supported", "single_source", "untested"),
#'                    c("B", "C", NA), c("FOXO3", "MIR23C", "NFE2L2"))
#' @export
#' @param family Logical vector from `consensus()`, TRUE where the label names
#'   a protein family or a complex rather than one gene. Such a label is a real
#'   answer, so it gets its own sentence: the collection names the group, not one
#'   gene, and no perturbation dataset can exist for it.
mechanism_writable <- function(verdict, level, regulator, sources = NULL,
                              family = NULL) {
    vapply(seq_along(verdict), function(i) {
        reg <- regulator[i]
        src <- if (is.null(sources)) character(0) else
            strsplit(sources[i] %||% "", ",", fixed = TRUE)[[1]]
        family_flag <- if (is.null(family)) logical(0) else family
        switch(
            verdict[i],
            supported = if (identical(level[i], "A") && any(src %in% c("tf", "mirna"))) {
                ## The strongest case the toolkit can make: a measured
                ## perturbation plus a collection that independently names the
                ## gene as a regulator of these targets.
                paste0(reg, ": a measured perturbation moves these genes and an ",
                       "independent target-set collection names it as their ",
                       "regulator. Report it as a candidate driver; the two ",
                       "evidence types agreeing is the strongest claim available.")
            } else if (identical(level[i], "A")) {
                paste0(reg, ": measured perturbation evidence supports a functional dependence.")
            } else {
                paste0(reg, ": target-set enrichment and at least one recorded mechanism agree; ",
                       "report the regulator as a candidate driver, verified experimentally.")
            },
            ## A single A-level source is still one source, but it is a
            ## measurement: saying "targets are enriched in target sets" about
            ## a knockdown would describe the wrong kind of evidence.
            single_source = if (isTRUE(family[i])) {
                ## A family name is a real answer, not a failure: the motif is
                ## enriched, but which member of the family drives it cannot be
                ## decided from target sets, and no knockdown exists for a
                ## complex. Saying so is more useful than a generic sentence.
                paste0(reg, ": target sets for this family or complex are ",
                       "enriched. The collection names the group, not one gene, ",
                       "so the member responsible is not identified here and no ",
                       "perturbation dataset can exist for it. Report the family ",
                       "and propose the member as the hypothesis.")
            } else if (mechanism_is_collision(reg)) {
                ## Say why, or a reader sees two sources listed and wonders why
                ## the verdict is not `supported`.
                paste0(reg, ": this label names different things in different ",
                       "sources (", mechanism_symbol_collisions()[[reg]], "), so ",
                       "their agreement is a coincidence of spelling. Do not ",
                       "combine them.")
            } else if (identical(level[i], "A")) {
                ## A knockdown moves the genes whether or not the perturbed
                ## gene regulates them: perturbing a downstream effector gives
                ## the same A grade as perturbing a driver. Only a collection
                ## that names the gene as a regulator distinguishes the two, so
                ## the sentence has to.
                if (any(src %in% c("tf", "mirna"))) {
                    paste0(reg, ": one measured perturbation moves these genes and ",
                           "its target sets are enriched too. A real dependence ",
                           "consistent with a regulatory role; corroborate before ",
                           "calling it a driver.")
                } else {
                    paste0(reg, ": one measured perturbation moves these genes, but ",
                           "nothing names it as a regulator -- it may be a ",
                           "downstream effector rather than an upstream driver. ",
                           "Do not report it as a mechanism.")
                }
            } else if (identical(sources_source(src), "mirna")) {
                ## A microRNA represses its targets, so the signature should be
                ## a down-regulation of genes carrying its seed. Without that
                ## direction the finding describes which genes carry the motif,
                ## not which microRNA is active.
                paste0(reg, ": predicted targets are over-represented, but a ",
                       "repressor's signature should appear among the repressed ",
                       "genes and this analysis cannot confirm that direction. ",
                       "Read it as a UTR-motif observation, not as microRNA activity.")
            } else {
                paste0(reg, ": targets are enriched in ", reg,
                       " target sets only. State it as a hypothesis, not as a mechanism.")
            },
            absent = paste0(reg, ": no supporting evidence found where evidence was available. ",
                            "An absence, not a refutation."),
            untested = paste0(reg, ": not examined. Do not write it as absent, and do not write it as supported."),
            NA_character_
        )
    }, "")
}


## A single label for the source of one writing decision.
sources_source <- function(src) {
    if (!length(src)) return(NA_character_)
    if (any(src == "mirna")) return("mirna")
    src[1]
}