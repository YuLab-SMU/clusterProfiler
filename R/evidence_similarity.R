#' Perturbation similarity from GPSAdb
#'
#' @description
#' Ask which measured perturbations produce a signature most like the caller's
#' gene list. This is the one question in the toolkit that does not start from a
#' prior: the target-set collections can only nominate factors someone already
#' annotated, while a similarity search compares the list against every
#' experiment GPSAdb holds.
#'
#' @details
#' The result is a ranking of *datasets*, not of regulators, and the distinction
#' is deliberate. "CEBPB knockdown looks like your list" is a statement about
#' an experiment; "CEBPB drives your list" is a statement about a mechanism, and
#' the first does not establish the second -- a knockdown of a downstream
#' effector produces a similar signature to a knockdown of its upstream
#' regulator. The ranking therefore lands in `$similar_perturbations` with the
#' datasets named, and does not enter `$regulators` where a verdict would be
#' read as a mechanism claim.
#'
#' What a similarity score can do is nominate. A factor that appears at the top
#' of this ranking and also has an enriched target set is a much better
#' candidate than either alone, and [consensus()] can see that once the caller
#' passes the names through `candidates`.
#'
#' @param genes Gene symbols.
#' @param log2fc Optional named numeric vector of fold changes, so the ranking
#'   can use direction rather than membership alone. Without it every gene is
#'   weighted equally and an opposite-direction experiment can rank highly.
#' @param limit Number of datasets to return.
#' @param min_genes GPSAdb refuses a ranked list of 200 or fewer genes; the
#'   server is the authority here, and its message is passed through.
#' @return A `data.frame` with `rank`, `dataset`, `perturbed_gene`, `method`,
#'   `system`, `score`, `effect` and `n_genes`.
#' @seealso [evidence()], [evidence_biotrigger()]
#' @examples
#' \dontrun{
#' evidence_similarity(c("HIF1A", "PGK1", "CA9"), limit = 10)
#' }
#' @export
evidence_similarity <- function(genes, log2fc = NULL, limit = 20,
                                min_genes = 200) {
    mechanism_require_mcp("gpsadb")
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(mechanism_empty_similarity())

    ## The server refuses short lists, and its own error explains why ("do not
    ## filter the list by p-value"). Checking here turns a round trip into a
    ## message that says what to do.
    if (length(genes) < min_genes) {
        rlang::abort(c(
            paste0("A similarity search needs more than ", min_genes,
                   " genes; this list has ", length(genes), "."),
            i = "GPSAdb rejects short lists because a handful of genes matches anything.",
            i = "Use the untrimmed result rather than a p-value-filtered subset."
        ))
    }

    ranked <- .similarity_input(genes, log2fc)
    res <- gpsadb::gpsadb_fastgpsa(ranked_list = ranked, limit = limit)
    d <- res$data
    if (is.null(d) || !nrow(d)) return(mechanism_empty_similarity())

    data.frame(
        rank        = as.integer(d$rank),
        dataset     = as.character(d$dataset_id),
        perturbed_gene = mechanism_alias_regulator(as.character(d$perturbed_gene)),
        method      = as.character(d$method %||% NA_character_),
        system      = as.character(d$cell_line %||% NA_character_),
        ## Higher score is a stronger match; the server's sign already encodes
        ## the direction of the effect.
        score       = as.numeric(d$gscore),
        effect      = tolower(as.character(d$effect %||% NA_character_)),
        n_genes     = as.integer(d$n_genes_in_dataset %||% NA_integer_),
        stringsAsFactors = FALSE
    )
}

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
    out[order(-abs(out$effect_size)), , drop = FALSE]
}

mechanism_empty_similarity <- function() {
    data.frame(rank = integer(), dataset = character(), perturbed_gene = character(),
               method = character(), system = character(), score = numeric(),
               effect = character(), n_genes = integer(), stringsAsFactors = FALSE)
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
.similarity_input <- function(genes, log2fc = NULL) {
    if (is.null(log2fc)) {
        return(data.frame(gene = genes, log2fc = seq_along(genes),
                          stringsAsFactors = FALSE))
    }
    if (is.null(names(log2fc))) {
        rlang::abort(c(
            "`log2fc` must be a named vector, with gene symbols as names.",
            i = "Without names it cannot be matched to `genes`."
        ))
    }
    v <- log2fc[genes]
    v[is.na(v)] <- 0
    data.frame(gene = genes, log2fc = as.numeric(v), stringsAsFactors = FALSE)
}
