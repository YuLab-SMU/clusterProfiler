#' Interpret an evidence bundle
#'
#' Hand an [evidence()] bundle to the interpretation pipeline, restricted to the
#' entities the evidence actually supports. This is the same call as
#' `interpret()` on an enrichment result -- `interpret()` dispatches on its
#' input -- so there is no second function to learn.
#'
#' The restriction is the point. Only entries that `adjudicate()` grades
#' `supported` or `single_source` **and** that carry level `A` or `B` reach the
#' prompt. A level C finding is a correlation: letting it into the prompt is how
#' a narrative ends up asserting a mechanism nobody demonstrated. Everything
#' excluded is still reported, in `$evidence_excluded`, so the omission is
#' visible rather than silent.
#'
#' `unsupported_claims` afterwards re-reads the model's prose and lists the
#' entities it named that were not in `$evidence_admitted`. It reuses
#' [evidence_status()] on each collection rather than on the bundle, because
#' `evidence_status()` expects an object with a `geneID` / `core_enrichment`
#' column and a bundle has neither.
#'
#' @param x An `evidenceBundle` from [evidence()].
#' @param ... Further arguments passed to the interpretation pipeline, such as
#'   `context`, `model` or `max_tokens`.
#' @return An `interpretation` object with two extra elements:
#'   `evidence_admitted` (the rows that reached the prompt) and
#'   `evidence_excluded` (the rows that did not, with the reason).
#' @seealso [evidence()], [adjudicate()], [interpret()]
#' @examples
#' \dontrun{
#' ev <- evidence(c("HIF1A", "PGK1", "CA9"))
#' interpret(ev)
#' }
#' @name interpret-evidence
NULL

## Turn a bundle into the named list of per-collection data frames that
## process_enrichment_input() already understands.
##
## The enrichment objects themselves are used when present, so the prompt keeps
## the same shape as a plain interpret() call; a bundle built with
## upstream = "tf" therefore looks exactly like interpret(list(tf = <enrichResult>)).
.process_evidence_input <- function(x, n_pathways = NULL, verdicts = NULL) {
    sources <- intersect(x$provenance$upstream, names(x))
    sources <- sources[vapply(sources, function(s) inherits(x[[s]], c("enrichResult",
                                                                     "gseaResult")), TRUE)]
    if (!length(sources)) return(NULL)

    admitted <- evidence_admitted(x, verdicts = verdicts)
    out <- vector("list", length(sources))
    names(out) <- sources

    for (s in sources) {
        obj <- x[[s]]
        keep <- admitted$term[admitted$source == s]
        df <- as.data.frame(obj)
        df <- df[!is.na(df$ID) & df$ID %in% keep, , drop = FALSE]
        if (nrow(df)) out[[s]] <- df else out[[s]] <- NULL
    }
    out <- out[!vapply(out, is.null, TRUE)]

    ## Measured evidence has no enrichment result to carry it: it is a table of
    ## (factor, dataset, direction). Dropping it would remove the only A-level
    ## evidence from the prompt, so a factor supported by a knockdown would be
    ## narrated as "predicted" by the target sets alone.
    meas <- mechanism_measured_prompt_rows(x, admitted)
    if (nrow(meas)) out[["perturbagen"]] <- meas

    if (!length(out)) return(NULL)
    out
}

## Measured evidence reshaped as the columns process_enrichment_input() reads.
mechanism_measured_prompt_rows <- function(x, admitted) {
    if (!"perturbagen" %in% names(x)) return(mechanism_empty_rows())
    keep <- admitted$term[admitted$source == "perturbagen"]
    if (!length(keep)) return(mechanism_empty_rows())

    m <- x$perturbagen
    m <- m[m$term %in% keep, , drop = FALSE]
    if (!is.null(m$padj_within_source) && nrow(m)) {
        m <- m[mechanism_significant(m), , drop = FALSE]
    }
    if (!nrow(m)) return(mechanism_empty_rows())

    ## The caller rbinds these with rows expanded from an enrichResult, so the
    ## column set has to match exactly; the columns that only make sense for an
    ## enrichment are filled with NA rather than omitted.
    data.frame(
        ID          = paste0(m$regulator, " (", m$direction, ")"),
        Description = paste0("perturbation of ", m$regulator, " [dataset ", m$term,
                             ", level A, direction ", m$direction, "]"),
        GeneRatio   = paste0(m$Count, "/", m$setSize),
        BgRatio     = NA_character_,
        RichFactor  = NA_real_,
        FoldEnrichment = NA_real_,
        oddsRatio   = NA_real_,
        zScore      = NA_real_,
        pvalue      = m$padj_within_source,
        p.adjust    = m$padj_within_source,
        qvalue      = NA_real_,
        geneID      = NA_character_,
        Count       = m$Count,
        stringsAsFactors = FALSE
    )
}

mechanism_empty_rows <- function() {
    mechanism_measured_prompt_rows_cols(0L)
}

mechanism_measured_prompt_rows_cols <- function(n) {
    data.frame(
        ID = character(n), Description = character(n), GeneRatio = character(n),
        BgRatio = character(n), RichFactor = numeric(n),
        FoldEnrichment = numeric(n), oddsRatio = numeric(n), zScore = numeric(n),
        pvalue = numeric(n), p.adjust = numeric(n), qvalue = numeric(n),
        geneID = character(n), Count = integer(n), stringsAsFactors = FALSE
    )
}

## The rows adjudicate() supports and the level allows.
evidence_admitted <- function(x, min_independent = 2, verdicts = NULL) {
    r <- x$regulators
    if (!nrow(r)) return(mechanism_empty_regulators())
    ## The caller's adjudication wins over a fresh one. Recomputing here with the
    ## default `min_independent` silently overrode a stricter decision the caller
    ## had already made -- which is the one thing adjudication as a separate step
    ## exists to prevent.
    cv <- .verdicts_or_adjudicate(x, min_independent, verdicts)

    keep <- cv$verdict %in% c("supported", "single_source") &
        !is.na(cv$max_level) & cv$max_level %in% c("A", "B")
    admitted_regulators <- cv$regulator[keep]

    r[r$regulator %in% admitted_regulators &
          r$named &
          mechanism_significant(r), , drop = FALSE]
}

## Abbreviations that a regulator narrative uses as words. Kept local to this
## check rather than added to .evidence_stop, which governs a different test.
.evidence_claims_stop <- c(
    ## Pathway and molecule abbreviations.
    "ECM", "GR", "ASM", "DEX", "GC", "IL", "TNF", "TGFB", "VEGF",
    "MAPK", "EGFR", "PI3K", "AKT", "MTOR", "COX", "NOS", "ROS",
    ## Cell lines and systems that appear in dataset descriptions. A reply
    ## quoting "knock_down CEBPB in BHK-21" is naming a cell line, not a gene.
    ## The list is long because GPSAdb spans hundreds of lines and any of them
    ## can surface in a reply about where the evidence came from.
    "BHK", "HELA", "MCF", "LNCAP", "A549", "HEK", "CHO", "COS", "K562",
    "HCT", "HEPG2", "SHRNA", "SIRNA", "KNOCKDOWN", "KNOCKOUT",
    "KYSE", "HUVEC", "MOLM", "BT", "EWS", "CCRF", "CEM", "ESF", "DU",
    "NPC", "HRMEC", "HCEC", "IPSC", "ASM", "HSAEC", "BEAS", "CALU", "H1299",
    "JURKAT", "THP1", "U937", "HL60", "NCI", "SW480", "HT29", "MDA", "SKOV",
    ## Units, methods and file words that a technical reply will use.
    "RNA", "DNA", "CDNA", "MRNA", "QPCR", "FDR", "BH", "GO", "KEGG",
    "GSEA", "TFT", "GMT", "MSIGDB", "STRING", "PPI")

mechanism_empty_regulators <- function() {
    data.frame(
        regulator = character(), source = character(), term = character(),
        named = logical(), level = character(), padj_within_source = numeric(),
        setSize = integer(), Count = integer(), family = logical(),
        evidence_kind = character(), direction = character(),
        dataset = character(), system = character(), method = character(),
        stringsAsFactors = FALSE
    )
}

## Everything that did not reach the prompt, with the reason.
evidence_excluded <- function(x, min_independent = 2, verdicts = NULL) {
    r <- x$regulators
    if (!nrow(r)) return(mechanism_empty_regulators()[0, ])
    cv <- .verdicts_or_adjudicate(x, min_independent, verdicts)
    ## The same table has to reach the admission decision too; recomputing it
    ## here would put the exclusion reason and the admitted set in disagreement.
    admitted <- evidence_admitted(x, min_independent = min_independent,
                                 verdicts = cv)
    left <- r[!r$term %in% admitted$term, , drop = FALSE]
    if (!nrow(left)) return(cbind(left, reason = character(0)))

    reason <- vapply(left$regulator, function(one) {
        v <- cv$verdict[cv$regulator == one]
        lv <- cv$max_level[cv$regulator == one]
        if (!length(v)) return("not graded")
        if (!v %in% c("supported", "single_source")) return(v)
        if (is.na(lv)) return("no level")
        paste0("level ", lv)
    }, "")
    cbind(left, reason = reason)
}

## Entities the narrative named that never reached the prompt and that the
## caller did not declare as background.
##
## The distinction matters: a model writing about a glucocorticoid experiment
## will name NR3C1, the receptor that defines the experiment, and that is context
## rather than an unevidenced claim. `design$background` names such entities so
## they are not reported; everything else the narrative introduces on its own is.
##
## Only well-formed symbols are reported. A regulator analysis invites prose
## like "AP-1" and "GR-NFAT", and a hyphenated bigram would otherwise be read as
## a gene symbol; abbreviations such as ASM or ECM are not symbols either. A
## token is reported only when it matches the shape of a HGNC symbol -- starts
## with a letter, continues with letters and digits, and is neither a
## hyphenated compound nor a bare acronym longer than five characters.
unsupported_claims <- function(interpretation, x, min_independent = 2,
                               verdicts = NULL) {
    ## The caller's adjudication is used for both steps. Computing it twice would
    ## let the admission decision and the exemptions below disagree about who was
    ## admitted, which is the same class of bug as re-adjudicating behind a
    ## caller's back.
    cv <- .verdicts_or_adjudicate(x, min_independent, verdicts)
    admitted <- evidence_admitted(x, min_independent = min_independent, verdicts = cv)
    background <- toupper(x$input$design$background %||% character(0))
    vocab <- unique(c(toupper(admitted$regulator), background,
                      toupper(unlist(strsplit(as.character(x$input$genes), "/")))))
    ## Being an input gene exempts a symbol from this check, because naming a
    ## gene that is in the data is not a claim about a regulator.
    ##
    ## That exemption is also how a rejected regulator walked through: the
    ## adjudication grades a symbol, declines to admit it, and the narrative
    ## names it as the mechanism -- and the check stays silent because the
    ## symbol is, separately, a differentially expressed gene. Measured on the
    ## MCF-7 example, where ESR1 is itself DE at lfc -1.94 and was knocked out.
    ##
    ## A symbol the adjudication graded and did not admit therefore keeps no
    ## exemption. Genes that were never graded are unaffected, so ordinary
    ## mentions of the input list stay out of the report.
    graded_out <- toupper(cv$regulator[!is.na(cv$verdict) &
        !cv$verdict %in% c("supported", "single_source")])
    vocab <- setdiff(vocab, graded_out)

    txt <- paste(unlist(interpretation[intersect(
        c("overview", "key_mechanisms", "regulatory_drivers", "crosstalk",
          "hypothesis", "narrative"), names(interpretation))]), collapse = " ")
    if (!nzchar(txt)) return(character(0))

    caps <- unique(unlist(regmatches(txt, gregexpr("\\b[A-Z][A-Z0-9]{1,11}\\b", txt))))
    ## No hyphens inside the token: 'AP-1' and 'GR-NFAT' are prose, not symbols.

    ## A family name for a regulator that WAS admitted is acceptable: writing
    ## STAT5 for an admitted STAT5A names the same family. The test is strictly
    ## one-directional, an admitted label must start with the token.
    ##
    ## The reverse direction (`startsWith(tok, vocab)`) must NOT be used: it
    ## accepted any real symbol whose prefix happens to be an admitted label,
    ## which made 1018 human symbols exempt, including STAT1, STAT3, FOXO1 and
    ## NFKB1 -- exactly the fabricated-regulator case this check exists to
    ## catch, so the check silently stopped working where it mattered most.
    ## Against the admitted labels only, not against `vocab`. `vocab` carries
    ## background and the input gene list, so a rejected regulator whose symbol
    ## prefixes an unrelated gene walked free: MYB was exempted because MYBL1 is
    ## an input gene, and MYBL1 is itself a rejected regulator in the MCF-7
    ## example. The exemption exists to let STAT5 stand for an admitted STAT5A,
    ## which is a statement about admitted labels.
    admitted_labels <- toupper(admitted$regulator)
    is_family_of_admitted <- function(tok) {
        any(startsWith(admitted_labels, tok))
    }
    ## A single letter followed by digits ('A4', 'Q4') is a list marker or a
    ## quality tier, never a symbol.
    caps <- caps[!grepl("^[A-Z][0-9]+$", caps)]
    unexplained <- setdiff(caps, c(vocab, .evidence_stop, .evidence_claims_stop))
    unexplained[!vapply(unexplained, is_family_of_admitted, TRUE)]
}

## The verdict table to use: the caller's if given, otherwise a fresh one.
##
## `verdicts = NULL` means "no adjudication was carried in", which is the case
## for the legacy entry points; those still get a verdict computed here rather
## than a different answer.
.verdicts_or_adjudicate <- function(x, min_independent, verdicts) {
    if (!is.null(verdicts)) return(verdicts)
    adjudicate(x, min_independent = min_independent)$verdicts
}
