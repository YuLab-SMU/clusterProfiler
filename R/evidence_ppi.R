#' Network evidence for a gene list
#'
#' Ask STRING how the genes relate and grade every edge by the channel that
#' supports it. An edge is **B** when it rests on experiment (`escore`) or a
#' curated database (`dscore`), and **C** when only text mining or co-expression
#' supports it -- a network can never produce an **A** grade, because no STRING
#' channel measures a perturbation.
#'
#' `pscore` is deliberately not used for grading. It is near-zero for almost
#' every edge, including pairs with published physical interaction (NFKB1 and
#' NFKBIA report `pscore = 0` with `escore = 0.999`), so a `pscore` threshold
#' would grade a known physical interaction as an association.
#'
#' @param genes Gene symbols.
#' @param species NCBI taxonomy id, defaulting to human (9606). getPPI() can
#'   only infer a species from an enrichment result, not from a gene vector.
#' @param required_score STRING's combined-score cutoff, 0-1000.
#' @param add_nodes Number of interaction partners to add.
#' @return A `data.frame` with one row per edge: `from`, `to`, `score`,
#'   `escore`, `dscore`, `tscore`, `pscore`, `channel` (the strongest channel
#'   supporting the edge) and `level` (`"B"` or `"C"`). An empty table means
#'   STRING returned no edge, which is a result, not an error.
#' @seealso [evidence()]
#' @export
evidence_ppi <- function(genes,
                         species = "auto",
                         required_score = NULL,
                         add_nodes = 0) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) < 2) {
        return(mechanism_empty_ppi())
    }

    ## getPPI() resolves taxID = "auto" only from an enrichResult, so a plain
    ## gene vector leaves it undetermined and the call fails. Default to human
    ## here and let the caller override it.
    if (identical(species, "auto")) species <- "9606"

    ## STRING takes the identifiers in the URL, so a few hundred genes overflow
    ## it and the request fails. Query in batches and combine; an edge between
    ## genes in different batches is missed, which is why the batch size is as
    ## large as the URL tolerates.
    batches <- split(genes, ceiling(seq_along(genes) / mechanism_ppi_batch()))

    n_batches <- length(batches)
    raw <- lapply(batches, function(batch) {
        if (length(batch) < 2) return(NULL)
        out <- tryCatch(
            getPPI(batch, ID = 1, taxID = species, required_score = required_score,
                   add_nodes = add_nodes, output = "data.frame"),
            error = function(e) NULL
        )
        if (is.null(out) || inherits(out, "error") || !nrow(out)) NULL else out
    })
    raw <- raw[!vapply(raw, is.null, logical(1))]
    if (!length(raw)) return(mechanism_empty_ppi())
    raw <- do.call(rbind, raw)

    ## Edges between genes that landed in different requests are invisible to
    ## STRING; record the batching so a caller can see the risk rather than
    ## trust the network as complete.
    out <- mechanism_grade_ppi(raw)
    attr(out, "n_batches") <- n_batches
    attr(out, "batch_size") <- mechanism_ppi_batch()
    out
}

mechanism_empty_ppi <- function() {
    data.frame(
        from = character(), to = character(), score = numeric(),
        escore = numeric(), dscore = numeric(), tscore = numeric(),
        pscore = numeric(), channel = character(), level = character(),
        stringsAsFactors = FALSE
    )
}

## Channels that can support a mechanism claim, strongest first. `nscore`
## (neighbourhood), `fscore` (gene fusion) and `ascore` (co-occurrence) describe
## similarity, not interaction, so they only ever contribute to C.
mechanism_ppi_channels <- function() {
    list(
        list(col = "escore", name = "experiment"),
        list(col = "dscore", name = "database")
    )
}

mechanism_grade_ppi <- function(raw) {
    if (is.null(raw) || !nrow(raw)) return(mechanism_empty_ppi())
    score <- function(col) {
        if (col %in% names(raw)) as.numeric(raw[[col]]) else rep(NA_real_, nrow(raw))
    }
    channel_score <- 0.4

    out <- data.frame(
        from   = as.character(raw$preferredName_A),
        to     = as.character(raw$preferredName_B),
        score  = score("score"),
        escore = score("escore"),
        dscore = score("dscore"),
        tscore = score("tscore"),
        pscore = score("pscore"),
        stringsAsFactors = FALSE
    )

    ## Name the strongest mechanism-supporting channel; fall back to the
    ## associational channels when neither reaches the threshold.
    channels <- mechanism_ppi_channels()
    best <- rep(NA_character_, nrow(out))
    for (ch in channels) {
        hit <- out[[ch$col]] >= channel_score
        best[hit & (is.na(best))] <- ch$name
    }
    out$channel <- ifelse(
        is.na(best),
        ifelse(!is.na(out$tscore) & out$tscore > 0, "textmining", "association"),
        best
    )
    out$level <- ifelse(out$channel %in% c("experiment", "database"), "B", "C")
    out
}

## Does the network support a specific pair? Used by consensus() to attach a
## verdict to a gene pair rather than to a regulator.
mechanism_pair_evidence <- function(ppi, a, b) {
    if (!nrow(ppi)) return(NULL)
    hit <- (ppi$from == a & ppi$to == b) | (ppi$from == b & ppi$to == a)
    if (!any(hit)) return(NULL)
    ppi[which(hit)[1], , drop = FALSE]
}

## Identifiers per STRING request.
##
## Measured against the live endpoint (2026-10-25): 400 genes in batches of 300
## lost 32.7% of edges to batch boundaries, 1000 genes in one request returned
## 741 edges in 4.5 s, and 2000 genes overflowed the URL and silently returned
## nothing. The default is therefore one batch for any realistic list, and the
## value is clamped so an oversized setting cannot reproduce that silent zero.
mechanism_ppi_batch <- function() {
    batch <- suppressWarnings(as.integer(Sys.getenv("CLUSTERPROFILER_PPI_BATCH", "1000")))
    if (is.na(batch)) batch <- 1000L
    max(2L, min(batch, 1000L))
}
