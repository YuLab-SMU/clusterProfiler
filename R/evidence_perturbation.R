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
#' @param system The biological system the datasets should come from, e.g.
#'   `"smooth muscle"` or `"HASM"`. Matched against each dataset's cell line and
#'   tissue, and used only to order the evidence: a dataset from the caller's
#'   system wins the per-regulator deduplication over a more significant one
#'   from elsewhere, and every row carries `in_system` so the match (or its
#'   absence) is visible rather than assumed. Set through `design$system`.
#' @param up,down The two halves of the caller's own contrast, as gene symbols.
#'   When given, `direction` is decided by a binomial test on the relevant half
#'   rather than by requiring every moved gene to go the same way -- a knockdown
#'   necessarily moves some genes up and some down, so that requirement was never
#'   met by real data and every dataset came back `mixed`. Restricted to one
#'   half, the question has an answer: among my up genes that this perturbation
#'   moves, do they move the way an activator's targets would? Both halves' counts
#'   and corrected p-values travel with the row, because which half is the right
#'   one to read depends on the direction the collection layer found, which is not
#'   known here. `NULL` (the default) keeps the older rule.
#' @param considered How many candidate regulators were eligible, when the
#'   caller already knows. Recorded on the result so a handful of datasets
#'   reads as "3 of the 30 strongest candidates had one" rather than "3
#'   factors drive this".
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
                                  considered = NULL,
                                  system = NULL,
                                  up = NULL,
                                  down = NULL) {
    mechanism_require_mcp("gpsadb")
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(mechanism_empty_perturbation())

    ## The sweep decides the answer, so say what it covered. With candidates
    ## the question is whether this regulator was perturbed, which is well
    ## posed; without them it is whether any sampled gene has a dataset, which
    ## is a property of the database rather than of the biology.
    ## A regulator is looked up under every name it is known by, because the
    ## name the collection uses is not always the name the database indexes.
    ## C3:TFT says `GR` where GPSAdb says `NR3C1`; querying only `GR` returns
    ## nothing, and a regulator with no dataset can never reach level A no
    ## matter what the literature says. See mechanism_gpsadb_names().
    query_genes <- if (is.null(candidates)) genes else
        unique(unlist(mechanism_gpsadb_names(as.character(candidates)),
                      use.names = FALSE))

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

    ## Which of these datasets were run in the caller's own system.
    datasets$in_system <- NA
    if (!is.null(system) && nzchar(system)) {
        datasets$in_system <- .gpsadb_in_system(datasets, system)
    }

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

        ## The direction is decided per half, by a test, not by requiring every
        ## sign to agree.
        ##
        ## The old rule asked whether all moved genes went the same way. A
        ## knockdown moves some genes up and some down, so on any list of a few
        ## hundred genes that is almost never true, and the answer was `mixed`
        ## for every dataset -- a property of the ruler, not of the data.
        ## Restricted to one half of the caller's own contrast, the question has
        ## an answer: among my up genes that this perturbation moves, do they
        ## move the way an activator's targets would?
        half_up <- .mechanism_sign_test(moved, up)
        half_down <- .mechanism_sign_test(moved, down)
        ## One correction across the two halves, because they are two tests of
        ## the same dataset and reporting either alone would overstate it.
        padj <- stats::p.adjust(c(half_up$p, half_down$p), method = "BH")
        half_up$padj <- padj[1]
        half_down$padj <- padj[2]
        half_up$direction <- .mechanism_half_direction(half_up)
        half_down$direction <- .mechanism_half_direction(half_down)

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
            ## Without a direction there is no half to test, so the old rule
            ## stands and says `mixed`: it is not wrong there, it is unaskable.
            direction          = if (!is.null(up) || !is.null(down)) half_up$direction else
                                  if (all(signs > 0)) "up" else
                                  if (all(signs < 0)) "down" else "mixed",
            up_rise            = half_up$rise,
            up_fall            = half_up$fall,
            up_sign_padj       = half_up$padj,
            down_rise          = half_down$rise,
            down_fall          = half_down$fall,
            down_sign_padj     = half_down$padj,
            dataset            = meta$title,
            system             = meta$system,
            method             = meta$method,
            ## Whether this dataset was run in the system the caller named.
            ## Carried rather than filtered on: a knockdown in BHK-21 still
            ## says the factor moves these genes, and dropping it would hide
            ## the caveat instead of stating it.
            in_system          = datasets$in_system[i],
            ## A perturbation dataset reports a p-value, so significance is
            ## decided by it rather than declared here.
            significant        = NA,
            effect_size        = NA_real_,
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
    ## One row per perturbed gene, because two datasets for the same gene are
    ## one claim, not two -- but which one to keep decides what the row means.
    ##
    ## Keeping the smallest p-value alone discarded the only airway smooth
    ## muscle dataset for SMAD3 (D20420, HCASMC) in favour of a prostate one
    ## whose p-value was smaller, so a caller who named their system never saw
    ## the dataset that matched it. The most relevant system therefore wins
    ## first, and significance breaks ties within it.
    ## The sort must be total, or the verdict is not reproducible. `padj` alone
    ## is not: a p-value underflows to exactly 0 for well-powered datasets, so
    ## SIX1's D25077 and D25078 tie, and `order()` then falls back on the order
    ## the rows happened to arrive in -- measured at 5 runs, the answer flipped
    ## between the two. The remaining keys settle it: more genes moved is the
    ## stronger claim, and the dataset id is unique, so no tie can survive to
    ## depend on row order.
    out <- out[order(-as.integer(out$in_system),
                     out$padj_within_source,
                     -out$Count,
                     out$dataset), , drop = FALSE]
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
        res$data[, intersect(c("dataset_id", "perturbed_gene", "cell_line", "tissue"),
                             colnames(res$data)), drop = FALSE]
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
        significant = logical(), effect_size = numeric(),
        in_system = logical(),
        ## The two halves' sign tests. Carried rather than collapsed: which half
        ## is the right one to read depends on the direction the collection
        ## layer found, and that is not known here.
        up_rise = integer(), up_fall = integer(), up_sign_padj = numeric(),
        down_rise = integer(), down_fall = integer(), down_sign_padj = numeric(),
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
## Only used where every slice is slow because of I/O. Serial is the default
## because the endpoint is not safe to query concurrently; see the note in the
## body. Forking is unavailable on Windows, which degrades to serial anyway.
mechanism_lapply <- function(x, fun) {
    ## Serial by default, and that is a correctness decision rather than a
    ## conservative one.
    ##
    ## The GPSAdb endpoint returns *different data* when it is queried
    ## concurrently, and says nothing about it: twenty parallel lookups of the
    ## same dataset produced the complete answer 16 times and a partial one 4
    ## times (231 of 379 significant genes, one run 53 of 91), while six serial
    ## lookups agreed exactly. Nothing distinguishes the partial answer -- same
    ## request, same status, no warning -- so a parallel sweep silently changes
    ## which dataset wins a regulator and therefore which verdicts come out.
    ## Measured end to end, the supported set flipped between runs.
    ##
    ## The sweep is I/O bound, so the parallel form is a speed-up that trades
    ## reproducibility for it. Raise `clusterProfiler_perturb_workers` to take
    ## that trade knowingly.
    workers <- getOption("clusterProfiler_perturb_workers", 1L)
    workers <- suppressWarnings(as.integer(workers))
    if (is.na(workers) || workers < 2L || length(x) < 2L ||
        .Platform$OS.type != "unix") {
        return(lapply(x, fun))
    }
    workers <- min(workers, length(x))
    ## Some hosts forbid forking -- R CMD check does, and so do sandboxes. The
    ## parallel form is a speed-up, so an environment that refuses it gets the
    ## serial answer rather than an error.
    ##
    ## mc.preschedule = FALSE spreads uneven requests rather than pinning a
    ## whole chunk to one child that happens to get the slow ones.
    tryCatch(
        parallel::mclapply(x, fun, mc.cores = workers, mc.preschedule = FALSE),
        error = function(e) lapply(x, fun)
    )
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

## Which datasets were run in the system the caller named.
##
## Matched on cell line and on tissue, case-insensitively and ignoring
## punctuation, because GPSAdb writes "HASMCs" where a caller would write "HASM"
## and "Smooth_Muscle" where they would write "smooth muscle". A caller naming a
## tissue therefore matches any of its cell lines, which is the useful
## direction: the tissue is what they know, the cell line is what the database
## records.
.gpsadb_in_system <- function(datasets, system) {
    norm <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
    want <- norm(system)
    if (!nzchar(want)) return(rep(NA, nrow(datasets)))

    haystack <- norm(paste(datasets$cell_line %||% "", datasets$tissue %||% ""))

    ## Both directions, element by element: the caller may name a tissue where
    ## the database records a cell line ("smooth muscle" against "HASMCs" is
    ## the case that matters), or a cell line where it records a cell line.
    ##
    ## vapply over the elements rather than grepl over the vector, because the
    ## second direction reverses pattern and subject and grepl would otherwise
    ## take only the first element as the pattern.
    vapply(haystack, function(h) {
        if (!nzchar(h)) return(FALSE)
        grepl(want, h, fixed = TRUE) || grepl(h, want, fixed = TRUE)
    }, logical(1), USE.NAMES = FALSE)
}

## One half's sign test.
##
## Restricted to the caller's own up (or down) genes, so the question is "among
## the genes I say went up, which way did this perturbation move them" rather
## than "did everything move together", which a knockdown never satisfies.
##
## Fewer than three genes gives no p-value. A binomial test on two observations
## can never reach 0.05, so it would be a test that always says `mixed` while
## looking like it had been run.
.mechanism_sign_test <- function(moved, half) {
    empty <- list(direction = NA_character_, rise = 0L, fall = 0L, p = NA_real_,
                  padj = NA_real_)
    if (is.null(half) || !length(half) || is.null(moved) || !nrow(moved)) {
        return(empty)
    }
    hit <- moved[moved$gene %in% half, , drop = FALSE]
    rise <- sum(hit$log2fc > 0)
    fall <- sum(hit$log2fc < 0)
    n <- rise + fall
    if (n < 3L) return(c(empty[c("direction", "p", "padj")],
                         list(rise = as.integer(rise), fall = as.integer(fall))))
    list(rise = as.integer(rise), fall = as.integer(fall),
         p = stats::binom.test(rise, n, 0.5)$p.value,
         padj = NA_real_, direction = NA_character_)
}

## Read a half's direction off its own test.
##
## `mixed` here means the test did not find a direction, which is not the same
## as the genes going both ways. Saying `mixed` for "no evidence of a side" is
## the honest report; saying `up` because one more gene rose would not be.
.mechanism_half_direction <- function(half, alpha = 0.05) {
    if (is.null(half$padj) || is.na(half$padj) || is.na(half$rise)) {
        return(NA_character_)
    }
    if (half$padj >= alpha) return("mixed")
    if (half$rise > half$fall) "up" else if (half$fall > half$rise) "down" else "mixed"
}
