## Narrative error injection: measure what the safety check actually catches.
##
## `interpret()` is bounded by the adjudication it was handed, in two places. The
## gate keeps unadmitted entities out of the prompt, and it has unit tests. The
## safety check -- `unsupported_claims()` -- runs on what the model wrote back,
## and had never been measured on its own.
##
## Measuring it by injecting into the narrative rather than the prompt is what
## makes the measurement exact: every error is one we placed, so the right answer
## is known, and no model call is needed. See dev_docs/mechanism-evaluation-spec.md
## for why only three of the six error kinds are in scope here.
##
## The three kinds this file covers:
##
##   C1  an unadmitted regulator named as the mechanism   -> must be caught
##   C3  a symbol that exists nowhere                     -> must be caught
##   C6  the untouched narrative                          -> must NOT be reported
##
## C2 (a single-source entity narrated as established), C4 (`untested` narrated
## as absent) and C5 (a direction stated backwards) are **not** measured, because
## `unsupported_claims()` is a symbol-vocabulary check and cannot see them. That
## they are missing from this report is a finding, not an oversight; the spec
## records them as gaps.

## The shape the safety check itself recognises: a letter followed by letters and
## digits, short enough to be a gene symbol. Candidates outside this shape (motif
## labels such as `AAANWWTGC_UNKNOWN`) cannot be caught by any symbol check, so
## scoring them would understate coverage for a reason that is not a defect.
.mechanism_symbol_shape <- function(x) {
    grepl("^[A-Z][A-Z0-9]{0,11}$", x)
}

## A narrative that is compliant by construction.
##
## Built from the verdict table so it names only admitted entities, which is what
## makes it usable as the C6 control: any report on this text is a false positive.
##
## Deliberately plain. It is a fixture, not a model of good writing -- a realistic
## model narrative would be more idiomatic and might contain words the check
## misreads, which would inflate the C6 false-positive rate for reasons that have
## nothing to do with the errors being injected.
.mechanism_narrative <- function(adj, max_regulators = 5L) {
    v <- adj$verdicts
    if (is.null(v) || !nrow(v)) return(NULL)
    keep <- v[v$verdict %in% c("supported", "single_source") & !is.na(v$regulator), ,
              drop = FALSE]
    if (!nrow(keep)) return(NULL)
    keep <- utils::head(keep[order(keep$padj), , drop = FALSE], max_regulators)

    each <- vapply(seq_len(nrow(keep)), function(i) {
        direction <- keep$direction[i]
        moved <- if (is.na(direction)) "differentially expressed" else
            paste0("", direction, "-regulated")
        sprintf("%s is %s here, and its target sets are enriched among the input genes.",
                keep$regulator[i], moved)
    }, "")
    paste(c("The following regulators are graded for this contrast.", each,
            "Single-source entries remain hypotheses rather than mechanisms."),
          collapse = " ")
}

## The regulators an injection may draw on: graded, not admitted, and shaped like
## a symbol so that the check has a fair chance.
##
## Regulators that are also input genes are kept in. Those are the ones that used
## to walk through untouched, and leaving them out would measure the check on the
## easy half of its own job.
.mechanism_injection_candidates <- function(adj) {
    v <- adj$verdicts
    if (is.null(v) || !nrow(v)) return(character(0))
    graded_out <- v$regulator[!is.na(v$verdict) &
        !v$verdict %in% c("supported", "single_source")]
    graded_out <- unique(graded_out[!is.na(graded_out)])
    graded_out <- graded_out[.mechanism_symbol_shape(toupper(graded_out))]
    setdiff(graded_out, c(.evidence_stop, .evidence_claims_stop))
}

## Place one error into a compliant narrative.
##
## Returns the corrupted text; the caller knows what it injected, which is the
## whole point. `C6` returns the narrative unchanged, so the same call path
## produces both the treatment and the control.
.mechanism_corrupt <- function(narrative, kind, symbol = NULL) {
    switch(kind,
        C1 = paste0(narrative, " ", symbol, " drives this programme directly."),
        C3 = paste0(narrative, " ", symbol, " acts upstream of the enriched terms."),
        C6 = narrative,
        stop("unknown injection kind: ", kind, call. = FALSE))
}

## What the safety check reports on a piece of text.
##
## A thin wrapper on purpose. This file measures `unsupported_claims()`; it does
## not reimplement it, or the measurement would be of something else.
.mechanism_detect <- function(narrative, adj) {
    if (is.null(adj) || is.null(adj$evidence)) return(character(0))
    unsupported_claims(list(narrative = narrative), adj$evidence,
                       verdicts = adj$verdicts)
}

## Run every in-scope kind and report per kind.
##
## Rates are reported separately and never combined. Merging them would hide the
## case this file exists to expose: a kind the check cannot see at all would
## disappear into an average that still looks respectable.
##
## `detected` counts cases where the injected symbol came back, not the total
## number of flagged symbols, so an unrelated flag cannot be scored as a catch.
.mechanism_injection_report <- function(adj, narrative = NULL, n_cases = 20L) {
    narrative <- narrative %||% .mechanism_narrative(adj)
    empty <- data.frame(kind = character(0), cases = integer(0),
                        injected = integer(0), detected = integer(0),
                        rate = numeric(0), stringsAsFactors = FALSE)
    if (is.null(narrative) || !nzchar(narrative)) return(empty)

    tally <- function(kind, symbols) {
        if (!length(symbols)) {
            return(data.frame(kind = kind, cases = 0L, injected = 0L,
                              detected = 0L, rate = NA_real_, stringsAsFactors = FALSE))
        }
        hits <- vapply(symbols, function(symbol) {
            text <- .mechanism_corrupt(narrative, kind, symbol)
            symbol %in% .mechanism_detect(text, adj)
        }, TRUE)
        data.frame(kind = kind, cases = length(symbols),
                   injected = length(symbols), detected = sum(hits),
                   rate = mean(hits), stringsAsFactors = FALSE)
    }

    ## C6 draws no symbol: the control is the narrative itself, and any report on
    ## it is a false positive.
    control <- data.frame(
        kind = "C6", cases = 1L, injected = 0L,
        detected = length(.mechanism_detect(narrative, adj)), rate = NA_real_,
        stringsAsFactors = FALSE)

    c1 <- utils::head(.mechanism_injection_candidates(adj), n_cases)
    ## Fabricated symbols are generated rather than drawn from data: nothing in
    ## the bundle can explain them, so a miss would be unambiguously the check's.
    ## Paired with C1 in count, so the two rates are comparable as the same amount
    ## of evidence rather than one being 20 cases and the other a hundred thousand.
    c3 <- paste0("ZQ", sprintf("%03d", seq_len(max(length(c1), 1L))))

    out <- rbind(tally("C1", c1), tally("C3", c3), control)
    out[order(match(out$kind, c("C1", "C3", "C6"))), , drop = FALSE]
}
