#' Decide which mechanisms the evidence actually supports
#'
#' Combine an `evidenceBundle` into one verdict per regulator, and return a
#' verdict object rather than a bare table. The point is not to rank by p-value
#' but to say which claims survive a change of evidence type, and -- just as
#' importantly -- which ones were never examined.
#'
#' **Adjudication is a step.** `evidence()` collects and grades, `adjudicate()`
#' decides, [interpret()] writes. Only the last one calls a model, and it needs
#' an `adjudication` -- not a bundle -- so the decision cannot be skipped on the
#' way to prose.
#'
#' The name is deliberate: this function does not find where sources agree, it
#' decides which claims survive. "Consensus" in bioinformatics means agreement
#' between clusterings, which is a different operation.
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
#'   regulator. Until then `adjudicate()` still does its job: it separates what
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
#' @return An S3 object of class `adjudication` with `verdicts`, `evidence`,
#'   `diagnostics`, `parameters` and `provenance`. `$verdicts` is a `data.frame`
#'   whose columns, in order, are `regulator`, `sources`, `n_sources`,
#'   `max_level`, `verdict`, `padj`, `family`, `evidence_kind`, `direction`,
#'   `measured_direction`, `perturbation_effect` and `writable` (the sentence
#'   the evidence permits; see [mechanism_writable()]).
#'
#'   `$provenance$evidence_sha` fingerprints the genes, the background, the
#'   collections and their release -- and deliberately not the time, so two runs
#'   over the same evidence hash alike and a comparison means "same evidence"
#'   rather than "same moment".
#'
#'   `direction` and `measured_direction` answer different questions and are
#'   kept apart deliberately. `direction` is the collection's: which way the
#'   regulator's programme moved in the caller's data, from the sign of the
#'   input half it was enriched in or the sign of its NES.
#'   `measured_direction` is the perturbation's: which way genes responded when
#'   the factor was knocked down, which is `mixed` for most real knockdowns
#'   because removing a factor moves some genes up and others down. Merging the
#'   two into one column made the knockdown's `mixed` overwrite the collection's
#'   established sign on exactly the best-supported entities.
#' @seealso [evidence()], [interpret()], [mechanism_writable()]
#' @examples
#' \dontrun{
#' ev <- evidence(setNames(c(2.1, 1.8, 1.4), c("PGK1", "CA9", "VEGFA")))
#' adj <- adjudicate(ev)
#' adj
#' consensus_findings(adj)
#' }
#' @export
adjudicate <- function(x, min_independent = 2, coverage_floor = NULL) {
    if (!inherits(x, "evidenceBundle")) {
        rlang::abort("`x` must be an evidenceBundle from evidence().")
    }
    if (is.null(coverage_floor)) coverage_floor <- x$input$design$coverage_floor

    r <- x$regulators
    if (!nrow(r)) {
        return(.adjudication(data.frame(
            regulator = character(), sources = character(), n_sources = integer(),
            max_level = character(), verdict = character(), padj = numeric(),
            family = logical(), evidence_kind = character(),
            direction = character(), measured_direction = character(),
            perturbation_effect = character(), writable = character(),
            stringsAsFactors = FALSE
        ), x, min_independent, coverage_floor))
    }

    ## A source that never examined these genes cannot report an absence.
    counts <- x$diagnostics$counts
    untested_source <- counts$source[
        !is.na(counts$coverage_floor) & counts$coverage_tested < counts$coverage_floor
    ]

    significant <- r[mechanism_significant(r), , drop = FALSE]
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
        ##
        ## But two sources can also be one witness. `perturbagen` and
        ## `biotrigger` are two views of the same GPSAdb experiment: one says
        ## "knocking this down changes these genes", the other "perturbing this
        ## matches this list". When both cite the same dataset, the agreement is
        ## the server agreeing with itself, and 17 of this package's `supported`
        ## verdicts were exactly that -- `C22orf46P` was level A supported on the
        ## strength of D25689 counted twice.
        datasets <- unique(mine$term[!is.na(mine$term) &
                                        mine$source %in% c("perturbagen", "biotrigger")])
        ## Distinct datasets, unless the sources name none (a collection term is
        ## not a dataset, so those do not collapse).
        n_witness <- if (length(datasets)) max(length(datasets), n) else n
        ## Overlapping dataset ids mean the sources are not independent even
        ## when they are different sources.
        shared_dataset <- length(datasets) < n &&
            all(mine$source %in% c("perturbagen", "biotrigger"))

        corroborated <- !shared_dataset && n_witness >= min_independent &&
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
            ## Two questions, two columns, because merging them destroyed both.
            ##
            ## `direction` is the collection's: which way the regulator's
            ## programme moved in the caller's data -- up-list enrichment or a
            ## positive NES. `measured_direction` is the perturbation's: which
            ## way the genes responded when the factor was knocked down. They
            ## are not the same measurement of the same thing, and taking the
            ## union of them meant a knockdown's `mixed` -- which is the normal
            ## answer, since removing a factor moves some genes up and others
            ## down -- overwrote the collection's direction on exactly the
            ## entities with the most evidence. Every one of the 8 regulators
            ## that only the GSEA path found came back `mixed`.
            direction  = .consensus_direction(mine$direction[
                !mine$source %in% c("perturbagen", "biotrigger")]),
            ## Reported on the same half `perturbation_effect` reads. Taking
            ## the perturbation row's own `direction` always reported the up
            ## half, so for a programme that runs down the two columns described
            ## different sets of genes.
            measured_direction = .consensus_measuring(mine)$measured,
            ## What the perturbation says the regulator DOES to the programme,
            ## when it says anything. A report, not a condition: requiring it
            ## would re-decide every verdict this package has already made, and
            ## those were not wrong -- they lacked one corroboration test, they
            ## did not rest on a false one.
            ##
            ## Deliberately not a logical "concordant" flag. On the hypoxia
            ## contrast a knockdown that pushes the up-genes further up reads as
            ## `FALSE` under that framing, which sounds like the evidence
            ## contradicts the finding; it does not, it says the factor represses
            ## them. Both are mechanisms, and collapsing them to true/false
            ## throws away the more interesting half.
            perturbation_effect = .consensus_perturbation_effect(mine),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out <- out[order(out$padj, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL
    out$writable <- mechanism_writable(out$verdict, out$max_level, out$regulator,
                                       sources = out$sources, family = out$family)
    .adjudication(out, x, min_independent, coverage_floor)
}

## The subset of a verdict table worth reading: entities the evidence actually
## speaks about. An `absent` verdict is reported once per source in
## `$diagnostics`; repeating it for every one of a thousand factor labels turns
## the sentence table into noise and hides the few that matter.
#' The findings an adjudication allows you to talk about
#'
#' Keep the entities whose verdict is `supported` or `single_source`, dropping
#' the ones that were examined and found wanting. This is the set a narrative
#' may discuss, and it is exported because the alternative -- filtering on the
#' verdict column by hand -- is easy to get subtly wrong: `absent` is evidence
#' against, `untested` is no evidence at all, and neither belongs in a claim.
#'
#' @param cv An `adjudication` from [adjudicate()], or its `$verdicts` table.
#' @return A `data.frame` of the admitted rows.
#' @seealso [adjudicate()], [collapse_findings()], [interpret()]
#' @examples
#' \dontrun{
#' consensus_findings(adjudicate(ev))
#' }
#' @export
consensus_findings <- function(cv) {
    ## Accepts the object as well as the frame: every caller has an
    ## adjudication, and reaching into $verdicts at each of them would spread
    ## knowledge of the object's shape across the package.
    cv <- adjudication_verdicts(cv)
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
#' @param verdict A `adjudicate()` verdict.
#' @param level An evidence level (`"A"`, `"B"`, `"C"`) or `NA`.
#' @param regulator The entity name.
#' @param sources Comma-separated source labels as `adjudicate()` records them.
#'   A measured perturbation is graded A whether it hit a regulator or a
#'   downstream effector, so the sentence for a single-source A finding depends
#'   on whether a target-set collection also names the gene as a regulator.
#' @return A character vector of sentences.
#' @seealso [adjudicate()]
#' @examples
#' mechanism_writable(c("supported", "single_source", "untested"),
#'                    c("B", "C", NA), c("FOXO3", "MIR23C", "NFE2L2"))
#' @export
#' @param family Logical vector from `adjudicate()`, TRUE where the label names
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

## Collapse a set of direction reports into one.
##
## `mixed` is not a value to average: it is what one source says when its own
## rows disagree, and a second `mixed` adds no information. A single consistent
## direction wins over an empty vector, and disagreement between sources is
## reported as disagreement rather than resolved silently.
.consensus_direction <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(NA_character_)
    if (length(unique(x)) == 1L) return(unique(x))
    "mixed"
}

## The collection programme, the perturbation's effect on it, and the method.
##
## One helper for both columns, because they must describe the same genes:
## `measured_direction` is what the perturbation did to the programme's own half,
## and `perturbation_effect` is how to read it. Computing them from different
## halves is how an activator came to be reported as a repressor on a knockdown
## the database reproduced exactly.
.consensus_measuring <- function(mine) {
    empty <- list(direction = NA_character_, measured = NA_character_,
                  sign = NA_character_)
    if (is.null(mine) || !nrow(mine)) return(empty)

    coll <- mine$direction[!mine$source %in% c("perturbagen", "biotrigger")]
    coll <- coll[!is.na(coll) & coll %in% c("up", "down")]
    if (!length(coll)) return(empty)
    direction <- unique(coll)
    ## Two collection rows pointing different ways means there is no single
    ## programme to compare against.
    if (length(direction) != 1L) return(empty)

    pert <- mine[mine$source %in% c("perturbagen", "biotrigger"), , drop = FALSE]
    if (!nrow(pert)) return(empty)
    sign <- mechanism_method_sign(pert$method[1])
    if (is.na(sign)) return(empty)

    ## Read the half the collection layer's direction refers to. The up half was
    ## tested into `direction` by the perturbation layer; the down half's verdict
    ## is derivable from the counts it stored.
    measured <- if (identical(direction, "up")) {
        pert$direction[1]
    } else {
        .mechanism_half_direction(list(rise = pert$down_rise[1],
                                       fall = pert$down_fall[1],
                                       padj = pert$down_sign_padj[1]))
    }
    if (is.null(measured) || is.na(measured) || !measured %in% c("up", "down")) {
        measured <- NA_character_
    }
    list(direction = direction, measured = measured, sign = sign)
}

## What the measured perturbation says the regulator does to the programme.
##
##   * the programme's genes fall when the factor is removed -> it promotes them
##   * they rise -> it restrains them
##
## "activator" and "repressor" are both mechanisms. Reporting a single
## concordant/discordant flag would have called the second one a failure, and on
## the hypoxia contrast that was 7 of 25 supported findings.
##
## NA when the comparison cannot be made -- no direction on one side, or the sign
## test found no side. "Could not check" and "checked, and it restrains them" are
## different answers.
.consensus_perturbation_effect <- function(mine) {
    m <- .consensus_measuring(mine)
    if (is.na(m$measured)) return(NA_character_)
    expected <- mechanism_expected_measured(m$sign)
    if (is.na(expected)) return(NA_character_)
    if (identical(m$measured, expected)) "activator" else "repressor"
}
## `called_at` and every other timestamp, so the same input hashes the same way
## on two runs and a SHA comparison means "same evidence", not "same moment".
.adjudication <- function(verdicts, bundle, min_independent, coverage_floor) {
    structure(
        list(
            verdicts    = verdicts,
            evidence    = bundle,
            diagnostics = bundle$diagnostics,
            ## What was asked, so a reader can reproduce the decision without
            ## reading the call that made it.
            parameters  = list(
                min_independent = min_independent,
                coverage_floor  = coverage_floor,
                species         = bundle$input$design$species,
                design          = bundle$input$design
            ),
            provenance  = list(
                evidence_sha = .adjudication_sha(bundle),
                call         = deparse(sys.call(-1L))[1],
                time         = as.character(Sys.time())
            )
        ),
        class = c("adjudication", "list")
    )
}

## A fingerprint of the evidence, not of the run.
##
## Everything in this list changes only when the input or the collections
## change. `bundles` are in a fixed order after the sort, so two runs over the
## same list hash identically however they were invoked.
.adjudication_sha <- function(bundle) {
    stable <- list(
        genes    = sort(bundle$input$genes),
        universe = sort(bundle$input$universe),
        ## Which half each gene is in, because that is what the directional
        ## layers read. The magnitudes are excluded on purpose: they decide no
        ## verdict, and including them would make the fingerprint change for a
        ## difference no answer can see.
        up       = sort(bundle$input$up %||% character(0)),
        down     = sort(bundle$input$down %||% character(0)),
        upstream = sort(bundle$provenance$upstream),
        release  = sort(bundle$provenance$release)
    )
    rlang::hash(stable)
}

#' Print an adjudication
#'
#' Shows what a reader needs in order to decide whether to go on to
#' [interpret()]: how many findings reached each verdict, and which regulators
#' the strongest ones name. The per-entity table is not printed -- at a few
#' hundred entities it buries the summary it is supposed to support.
#'
#' @param x An `adjudication` from [adjudicate()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.adjudication <- function(x, ...) {
    v <- x$verdicts
    cat("<adjudication>\n")
    cat("  evidence:     ", length(x$evidence$input$genes), " genes, ",
        paste(x$evidence$provenance$upstream, collapse = ", "), "\n", sep = "")

    counts <- table(factor(v$verdict, levels = c("supported", "single_source",
                                                 "absent", "untested")))
    cat("  verdicts:     ", paste0(names(counts), " ", as.integer(counts),
                                  collapse = "  "), "\n", sep = "")

    ## A source that could not be consulted, or that did not cover enough of the
    ## list, changes what `untested` means -- say so here rather than leaving it
    ## to be discovered in $diagnostics.
    floors <- vapply(x$diagnostics$counts$status %||% character(0),
                     function(s) identical(s, "untested"), TRUE)
    if (any(floors)) {
        cat("  coverage:      ", sum(floors), " source(s) below the coverage floor",
            " -- their entities are untested, not absent\n", sep = "")
    }

    sup <- v[v$verdict == "supported", , drop = FALSE]
    if (nrow(sup)) {
        sup <- utils::head(sup, 3)
        cat("  supported:    ",
            paste0(sup$regulator, " (", sup$max_level, ")",
                   ifelse(is.na(sup$perturbation_effect), "",
                          paste0(", ", sup$perturbation_effect)),
                   collapse = "; "), "\n", sep = "")
    }
    ## A supported finding resting on a few percent of the input is not a weaker
    ## claim by rule, but it is a different claim to read. Named explicitly
    ## rather than folded into the verdict, which stays as it was.
    thin <- .mechanism_thin_supported(x)
    if (!is.null(thin)) {
        cat("  thin:         ", paste(thin, collapse = "; "),
            " -- supported, but the measured evidence covers little of the input\n", sep = "")
    }
    cat("  sha:          ", substr(x$provenance$evidence_sha, 1, 16), "\n", sep = "")
    invisible(x)
}

## The verdict table, for callers that only want the frame.
##
## `consensus_findings()` and several helpers took the data.frame directly, and
## making them reach into `$verdicts` at every call site would spread knowledge
## of the object's shape across the package.
adjudication_verdicts <- function(x) {
    if (inherits(x, "adjudication")) x$verdicts else x
}
