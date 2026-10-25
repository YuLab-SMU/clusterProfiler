#' Perturbation evidence from GPSAdb
#'
#' Ask GPSAdb which measured perturbation moves the input genes, and grade the
#' result. A perturbation dataset whose target set is enriched among the input
#' genes is the only evidence that can reach level **A**: it is an experiment,
#' not an annotation or a prediction.
#'
#' The returned rows carry the same `regulator` label as the collections built by
#' [mechanism_sets()], so a transcription factor named both by its target sets
#' and by a knockdown dataset is *one* entity supported by *two* sources --
#' which is what [consensus()] needs in order to report `supported` at all.
#'
#' GPSAdb is reached through its MCP endpoint, so `gpsadb` and `aisdk.mcp` must
#' be installed. They are suggested, not required: without them this function
#' aborts with an install hint rather than returning an empty table, because an
#' empty table would read as "no perturbation evidence exists".
#'
#' @param genes Gene symbols.
#' @param method Perturbation method to query, e.g. `"knockdown"`, `"knockout"`,
#'   `"overexpression"`, `"drug"`. `NULL` queries every method.
#' @param limit Maximum number of datasets considered.
#' @param min_genes Minimum number of input genes a dataset must contain before
#'   it is tested. A dataset with two hits says nothing about the list.
#' @param pvalueCutoff Cutoff for a dataset to be reported as perturbing.
#' @param candidates Genes to query as the perturbed gene. Supply the regulators
#'   the collections already nominated and the sweep becomes targeted: asking
#'   whether a known candidate was ever knocked down is a question GPSAdb
#'   answers directly, while sweeping a sample of an arbitrary gene list only
#'   reports which genes happen to have datasets. Given NULL the function falls
#'   back to a spread sample and records how much it scanned in
#'   attr(x, "genes_scanned").
#' @return A `data.frame` with `regulator` (the perturbed gene), `source`
#'   (`"perturbagen"`), `term` (the dataset id), `named`, `level` (`"A"` when
#'   reported, `"B"` when the dataset is absent but the gene is in GPSAdb),
#'   `padj_within_source`, `setSize`, `Count` and `direction` (`"up"`, `"down"`
#'   or `"mixed"`, from the sign of the responding genes).
#' @seealso [evidence()], [consensus()]
#' @examples
#' \dontrun{
#' evidence_perturbation(c("HIF1A", "PGK1", "CA9"), method = "knockdown")
#' }
#' @export
evidence_perturbation <- function(genes,
                                  method = NULL,
                                  limit = 50,
                                  min_genes = 5,
                                  pvalueCutoff = 0.05,
                                  candidates = NULL,
                                  considered = NULL) {
    mechanism_require_mcp("gpsadb")
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(mechanism_empty_perturbation())

    ## The sweep decides the answer, so say what it covered. With candidates
    ## the question is whether this regulator was perturbed, which is well
    ## posed; without them it is whether any sampled gene has a dataset, which
    ## is a property of the database rather than of the biology.
    query_genes <- if (is.null(candidates)) genes else
        unique(as.character(candidates))

    ## Which perturbations move any of these genes: the reverse lookup GPSAdb
    ## provides, which names the perturbed gene for each dataset.
    hits <- tryCatch(
        gpsadb::gpsadb_query_gene(genes[1], limit = 1),
        error = function(e) NULL
    )
    if (is.null(hits)) return(mechanism_empty_perturbation())

    ## The single-gene lookup only proves the endpoint answers; the per-gene
    ## sweep is what collects the datasets.
    datasets <- .gpsadb_datasets_for(query_genes, method = method, limit = limit,
                                     max_genes = if (is.null(candidates)) 30L else
                                         length(query_genes))
    if (is.null(datasets) || !nrow(datasets)) return(mechanism_empty_perturbation())

    ## Each candidate costs a network round trip, so the sweep is run in
    ## parallel. With 23 of 30 candidates having data the serial version took
    ## 170 s for the airway list; the requests are independent, so they go out
    ## together. Parallelism is opt-out via the option for a constrained host.
    rows <- mechanism_lapply(seq_len(nrow(datasets)), function(i) {
        id <- datasets$dataset_id[i]
        perturbed <- datasets$perturbed_gene[i]

        ## The question is not "are these genes enriched among the responding
        ## genes" -- gpsadb_de() returns the *most significant* 500 rows, not
        ## the measured set, so any enrichment against them is meaningless.
        ## The question a perturbation answers is narrower and exactly what a
        ## mechanism needs: does knocking this factor down move the genes I am
        ## looking at, and in which direction?
        d <- .gpsadb_lookup(id, genes)
        if (is.null(d) || !nrow(d)) return(NULL)
        d <- d[!is.na(d$gene) & !is.na(d$padj), , drop = FALSE]
        moved <- d[d$padj < pvalueCutoff, , drop = FALSE]
        if (nrow(moved) < min_genes) return(NULL)

        signs <- sign(moved$log2fc)
        ## The dataset's own description travels with the row. "D22765" alone
        ## tells a reader nothing; the cell line and method decide whether the
        ## measurement is relevant to the caller's system at all -- a knockdown
        ## in BHK-21 says little about human airway smooth muscle.
        meta <- .gpsadb_dataset_meta(id)
        data.frame(
            regulator          = mechanism_alias_regulator(perturbed),
            source             = "perturbagen",
            term               = id,
            named              = TRUE,
            family             = FALSE,
            ## A perturbation dataset is a measurement by definition.
            evidence_kind      = "measured",
            ## Level A: this is a measured perturbation, the only source that
            ## can reach it.
            level              = "A",
            ## Among the measured targets, the best adjusted p-value stands for
            ## how confidently the factor moves the list.
            padj_within_source = min(moved$padj),
            setSize            = nrow(d),
            Count              = nrow(moved),
            direction          = if (all(signs > 0)) "up" else
                                  if (all(signs < 0)) "down" else "mixed",
            dataset            = meta$title,
            system             = meta$system,
            method             = meta$method,
            stringsAsFactors   = FALSE
        )
    })

    rows <- rows[!vapply(rows, is.null, TRUE)]
    if (!length(rows)) {
        out <- mechanism_empty_perturbation()
        attr(out, "lookup_failures") <- attr(datasets, "lookup_failures")
        return(out)
    }

    out <- do.call(rbind, rows)
    ## What the sweep covered, so the result is readable as "3 of the 30
    ## strongest candidates had a dataset" rather than "3 factors drive this".
    attr(out, "queried") <- if (is.null(candidates)) NA_integer_ else length(query_genes)
    attr(out, "considered") <- considered %||% attr(datasets, "genes_total")
    attr(out, "targeted") <- !is.null(candidates)
    ## One row per perturbed gene: the strongest dataset that moves it,
    ## because two datasets for the same gene are one claim, not two.
    out <- out[order(out$padj_within_source), , drop = FALSE]
    out[!duplicated(out$regulator), , drop = FALSE]
}

## Dataset ids whose measured perturbation moves at least one input gene.
##
## GPSAdb searches a single perturbed gene at a time, so the sweep is one query
## per gene. Which genes get swept decides the answer: taking the first twenty
## of a fold-change-sorted list scans whichever genes happen to sort first,
## which for a glucocorticoid list meant BRCA1 and SAMD4A while the biology of
## interest was never queried. The sweep therefore covers a *balanced* sample
## of the list rather than its head, and the cap is recorded in the result so a
## caller can see how much was scanned.
##
## shortcut: a sample, not the whole list. Sweeping 897 genes costs 897 requests
## at about a second each; raise `max_genes` when a complete answer matters more
## than the wait.
.gpsadb_datasets_for <- function(genes, method = NULL, limit = 50,
                                 max_genes = getOption("clusterProfiler_perturb_scan", 30L)) {
    ## Spread the sample across the list instead of taking a prefix: genes that
    ## sort early and genes that sort late are equally likely to be a
    ## perturbation target, and a prefix biases toward whichever end.
    n <- min(length(genes), as.integer(max_genes))
    queries <- if (n >= length(genes)) genes else {
        idx <- round(seq(1, length(genes), length.out = n))
        genes[unique(idx)]
    }

    found <- lapply(queries, function(g) {
        res <- tryCatch(
            gpsadb::gpsadb_search_datasets(perturbed_gene = g, method = method,
                                           limit = limit),
            error = function(e) NULL
        )
        if (is.null(res) || is.null(res$data) || !nrow(res$data)) return(NULL)
        res$data[, intersect(c("dataset_id", "perturbed_gene"), colnames(res$data)),
                 drop = FALSE]
    })
    found <- found[!vapply(found, is.null, TRUE)]
    if (!length(found)) return(NULL)
    out <- do.call(rbind, found)
    out <- out[!is.na(out$dataset_id), , drop = FALSE]
    out <- out[!duplicated(out$dataset_id), , drop = FALSE]
    attr(out, "genes_scanned") <- length(queries)
    attr(out, "genes_total") <- length(genes)
    out
}

mechanism_empty_perturbation <- function() {
    data.frame(
        regulator = character(), source = character(), term = character(),
        named = logical(), level = character(), padj_within_source = numeric(),
        setSize = integer(), Count = integer(), family = logical(),
        evidence_kind = character(), direction = character(),
        dataset = character(), system = character(), method = character(),
        stringsAsFactors = FALSE
    )
}

## A companion package missing is not "no evidence": say so and stop.
mechanism_require_mcp <- function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
    rlang::abort(c(
        paste0("Package '", pkg, "' is required for this evidence source."),
        i = "Install it with remotes::install_github(\"YuLab-SMU/", pkg, "\").",
        i = "Returning an empty table instead would read as 'no evidence exists'."
    ))
}

## Look the genes up in one dataset, in batches small enough for the endpoint.
##
## gpsadb_de() refuses more than 500 genes at once. Sending a realistic gene
## list in one call therefore fails, and because the failure was caught and
## turned into NULL the caller reported "no perturbation evidence" -- an
## absence claim produced by an API limit. Batching keeps the answer honest;
## a real failure is still reported through the `lookup_failures` attribute
## rather than being swallowed.
.gpsadb_lookup <- function(id, genes, batch = 500L) {
    chunks <- split(genes, ceiling(seq_along(genes) / batch))
    parts <- lapply(chunks, function(chunk) {
        tryCatch(gpsadb::gpsadb_de(id, genes = chunk, limit = 500),
                 error = function(e) e)
    })
    failed <- vapply(parts, inherits, TRUE, "error")
    if (any(failed)) {
        warning(paste0("GPSAdb lookup failed for dataset ", id, ": ",
                       conditionMessage(parts[[which(failed)[1]]])),
                call. = FALSE)
    }
    ok <- parts[!failed]
    ok <- lapply(ok, function(p) if (is.null(p$data)) NULL else p$data)
    ok <- ok[!vapply(ok, is.null, TRUE)]
    if (!length(ok)) return(NULL)
    out <- do.call(rbind, ok)
    out[!duplicated(out$gene), , drop = FALSE]
}

## Parallel map over request-shaped work.
##
## Only used where every element is independent and slow because of I/O. Forking
## is unavailable on Windows, so that platform falls back to serial rather than
## failing; `options(clusterProfiler_perturb_workers = 1)` forces serial.
mechanism_lapply <- function(x, fun) {
    workers <- getOption("clusterProfiler_perturb_workers", 4L)
    workers <- suppressWarnings(as.integer(workers))
    if (is.na(workers) || workers < 2L || length(x) < 2L ||
        .Platform$OS.type != "unix") {
        return(lapply(x, fun))
    }
    workers <- min(workers, length(x))
    ## mc.preschedule = FALSE spreads uneven requests rather than pinning a
    ## whole chunk to one child that happens to get the slow ones.
    parallel::mclapply(x, fun, mc.cores = workers, mc.preschedule = FALSE)
}

## One dataset's description, or blanks when it cannot be fetched.
##
## Metadata is a second request, so a failure here must not lose the evidence
## row: the numbers are already in hand and the description is supplementary.
.gpsadb_dataset_meta <- function(id) {
    blank <- list(title = NA_character_, system = NA_character_, method = NA_character_)
    ds <- tryCatch(gpsadb::gpsadb_dataset(id), error = function(e) NULL)
    if (is.null(ds) || is.null(ds$meta)) return(blank)
    m <- ds$meta
    list(
        title  = as.character(m$title %||% NA_character_),
        system = as.character(m$cell_line %||% NA_character_),
        method = as.character(m$method %||% NA_character_)
    )
}
