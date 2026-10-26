#' Upstream regulators from GPSAdb's biotrigger
#'
#' @description
#' Ask which perturbed genes move the caller's list, which is the reverse
#' question to [evidence_similarity()]: instead of "which experiment looks like
#' this list" it asks "perturbing what changes this list".
#'
#' @details
#' A biotrigger hit is a measured perturbation whose effect overlaps the list,
#' so it is graded **B**, the same as a GTRD target set: the factor was
#' measured, but the measurement is of a signature, not of binding, and the
#' overlap of two signatures is weaker than a direct assay of the same factor.
#' It is not graded C, because unlike a motif it is not a prediction.
#'
#' The genes returned are *perturbed genes*, not necessarily regulators. A hit
#' on a metabolic enzyme means perturbing that enzyme shifts the list, which is
#' a downstream relationship. The returned rows carry the effect direction so
#' the caller can judge, and a hit is not a claim that the gene regulates
#' anything.
#'
#' @param genes Gene symbols.
#' @param effect Restrict to `"activated"`, `"inhibited"` or `"both"`.
#' @param limit Number of hits.
#' @param min_genes Minimum list size worth asking about.
#' @return A `data.frame` with `regulator`, `source`, `term`, `named`, `level`,
#'   `padj_within_source`, `setSize`, `Count`, `family`, `evidence_kind`,
#'   `direction` and the dataset columns, in the shape [evidence()] uses.
#' @seealso [evidence_similarity()], [evidence()]
#' @examples
#' \dontrun{
#' evidence_biotrigger(c("HIF1A", "PGK1", "CA9", "BNIP3", "LDHA"))
#' }
#' @export
evidence_biotrigger <- function(genes, effect = c("both", "activated", "inhibited"),
                                limit = 20, min_genes = 5) {
    mechanism_require_mcp("gpsadb")
    effect <- match.arg(effect)
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) < min_genes) return(mechanism_empty_biotrigger())

    res <- gpsadb::gpsadb_biotrigger(genes, effect = effect, limit = limit)
    d <- res$data
    if (is.null(d) || !nrow(d)) return(mechanism_empty_biotrigger())

    out <- data.frame(
        regulator  = mechanism_alias_regulator(as.character(d$perturbed_gene)),
        source     = "biotrigger",
        term       = as.character(d$dataset_id),
        named      = TRUE,
        ## Measured, but weaker than a direct perturbation assay: the evidence
        ## is the overlap of two signatures.
        level      = "B",
        ## GPSAdb reports an effect size, not a p-value. The magnitude of the
        ## enrichment score stands in for strength, and a p-value would be a
        ## number the server never produced.
        padj_within_source = NA_real_,
        setSize    = as.integer(d$n_genes_in_dataset %||% NA_integer_),
        Count      = as.integer(d$n_genes_in_dataset %||% NA_integer_),
        family     = mechanism_is_family_label(
            mechanism_alias_regulator(as.character(d$perturbed_gene))),
        evidence_kind = "measured",
        direction  = ifelse(tolower(as.character(d$effect)) == "activated",
                            "up", "down"),
        dataset    = as.character(d$dataset_id),
        system     = as.character(d$cell_line %||% NA_character_),
        method     = as.character(d$method %||% NA_character_),
        stringsAsFactors = FALSE
    )
    ## The enrichment score decides ordering; it is carried through so the
    ## caller can see why one hit ranked above another.
    out$effect_size <- as.numeric(d$es)
    ## The server reports a score, not a test. A hit is significant by the
    ## server's own construction -- it returned it as a hit -- so the row says
    ## so itself rather than leaving consensus() to read a missing p-value as
    ## "not significant", which is what a padj-only rule would do.
    out$significant <- TRUE
    out[order(-abs(out$effect_size)), , drop = FALSE]
}



mechanism_empty_biotrigger <- function() {
    out <- mechanism_empty_perturbation()
    out$effect_size <- numeric()
    out
}

## The ranked list GPSAdb wants: one row per gene, in the order the caller gave.
##
## Direction matters more than membership -- an experiment that moves the same
## genes the other way is the opposite finding -- so a fold change is used when
## one is supplied and the caller's order otherwise.
