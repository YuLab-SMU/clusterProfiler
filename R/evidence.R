#' Grade upstream evidence for a gene list
#'
#' Turn a gene list into a table of candidate mechanisms, each carrying an
#' evidence level that says how strongly it may be stated. This is the
#' mechanism-exit counterpart of a pathway enrichment: the collections tested
#' name *upstream participants* (transcription factors, microRNAs, ...) instead
#' of biological processes, so an enriched set answers "what could have
#' produced these genes" rather than "what do these genes do".
#'
#' Each requested collection is tested separately against the same gene list,
#' and each is corrected within itself. Pooling them would let a large
#' collection dilute a small one: C3:TFT holds about 1100 sets and C3:MIR about
#' 2600, with different target definitions, so their p-values are not on one
#' scale. `$regulators` therefore carries `padj_within_source`, never a single
#' shared `padj`.
#'
#' Every returned collection must be graded the same way: an over-represented
#' target set is evidence that the regulator's *fingerprint* sits in the list,
#' not evidence that the regulator caused it. `$regulators$level` says how far
#' a claim may go, and `$diagnostics` says how much of the list the collection
#' could speak about at all.
#'
#' @param genes Gene symbols, or a **named numeric vector** of gene symbols and
#'   their statistics. A named numeric is tested by direction: the positive half
#'   and the negative half are run separately against each collection, so a
#'   regulator's row carries whether its programme rose or fell. Measured on the
#'   hypoxia contrast this took `HIF1_Q3` from padj 1.2e-4 (both directions
#'   mixed into one list) to 1.0e-5, because asking "which set is enriched among
#'   what changed" dilutes the signal when half the changes go the other way.
#'   Both halves are corrected together, since splitting doubles the number of
#'   hypotheses. An `enrichResult` (its
#'   `geneID` column is expanded), a `gseaResult` (the union of its
#'   `core_enrichment` genes, with the ranked list as the background), or a
#'   `data.frame` with a column named by `design$gene_col`.
#'
#'   A `gseaResult` contributes the genes that produced each enriched set's
#'   ranking rather than every gene in `@geneList`: asking which regulator
#'   explains all of a transcriptome has no answer, and on the airway list that
#'   form returned no `supported` regulator at all. Measured against known
#'   answers, the GSEA form and the `enrichResult` form rank the true regulator
#'   identically (MRR 0.046 against 0.047), because the limit is the
#'   over-representation test's power at small overlaps, not the input.
#' @param pathway The conventional enrichment layer. `FALSE` (the default)
#'   skips it. `TRUE` runs it here, on the up and the down half separately, and
#'   then asks which regulator's target set explains each pathway -- a different
#'   question from the collection layer, which asks which regulator's target set
#'   is over-represented in the whole list. An `enrichResult` or `gseaResult`
#'   reuses a result the caller already built instead of recomputing it, which is
#'   also the way to pass one: `evidence()` no longer takes a result as `genes`,
#'   because reading genes back out of it silently drops every gene that appears
#'   only in a pathway below its p-value cutoff.
#'
#'   The pathway rows are kept in `$pathways`, never in `$regulators`. They are
#'   derived from the same collections, so counting them as a source would let
#'   one collection corroborate itself.
#' @param design A list describing the experiment; see Details.
#' @param upstream Collections to test: `"tf"`, `"tfunion"`, `"mirna"`,
#'   `"perturbagen"`, `"cellstate"` or `"tissue"`. Defaults to `"tf"` alone.
#'
#'   `"tfunion"` is C3:TFT plus four public transcription-factor libraries
#'   (ENCODE, ChEA, TF_Perturbations and TRRUST), which together name 315
#'   factors C3:TFT has no set for -- including `EPAS1`, `NFKB1`, `TP53`,
#'   `CTCF`, `ESR1`, `JUN` and `FOS`. Prefer it when the question is "which
#'   regulator drives this list" rather than "which of these specific motif
#'   families is enriched": on a hypoxia contrast it recovered the paper's
#'   answer (HIF-2) where `"tf"` could not, at 2869 sets against 1114, and
#'   three random gene sets of the same size still returned no significant
#'   set. See [mechanism_tfunion()].
#'
#'   Pass it **together with a measured source**, e.g.
#'   `c("tfunion", "perturbagen")`. The union is one source however many sets
#'   it holds, and corroboration counts sources rather than rows, so on its own
#'   it can only reach `single_source` -- a page of those reads like a failure
#'   when the finding is simply uncorroborated.
#'
#'   `"mirna"` is
#'   off by default for a measured reason: its sets are predicted seed matches,
#'   so an enrichment describes which genes carry the motif rather than which
#'   microRNA is active. Measured on the airway list, 7 miRNAs were
#'   significant among the UP-regulated genes and 0 among the down ones -- the
#'   opposite of what a repressor predicts -- so the direction that would make
#'   the finding interpretable is absent. Ask for it explicitly when the
#'   question is about UTR composition.
#' @param universe Background genes. `NULL` uses every symbol in `org.Hs.eg.db`.
#' @param tft Which part of the transcription-factor collection to test. `NULL`
#'   (the default) uses all of C3:TFT. `"gtrd"` keeps only the ChIP-defined sets,
#'   whose identifiers are gene symbols: measured on the airway list that takes
#'   the share of labels naming a family or a consensus sequence from 40% to 1%,
#'   and those labels are the ones no perturbation dataset can ever answer for.
#'   `"motif"` selects the predicted half instead; a character vector selects
#'   several, e.g. `c("gtrd", "motif")`.
#'
#'   Only `upstream = "tf"` honours this. Asking for it together with
#'   `"tfunion"` is an error rather than a silent no-op: the union's other four
#'   libraries have no GTRD/motif split, so ignoring the filter would answer a
#'   wider question than the one asked.
#' @param network Whether to query STRING for interaction evidence. Off by
#'   default because it needs the network; see [evidence_ppi()] for the grading.
#' @param minGSSize Minimum gene-set size. The default is 5, not `enricher()`'s
#'   10, so that small motif sets are tested; `maxGSSize` is `Inf` because the
#'   largest GTRD target sets carry more than a thousand genes and the default
#'   upper bound of 500 would drop exactly the broad regulators that matter.
#' @param maxGSSize Maximum gene-set size.
#' @param print_pvalue,print_qvalue Cutoffs used by the print method only. They
#'   do not filter the result: `$regulators` always holds every tested set, so
#'   that a set filtered out for display can still be distinguished from a set
#'   that was never tested.
#' @return An `evidenceBundle`: a list with `input`, `sets`, one element per
#'   tested collection (`tf`, `mirna`), `regulators`, `xref`, `diagnostics` and
#'   `provenance`.
#'
#'   `diagnostics` reports counts whose denominators differ and must not be
#'   conflated: `n_input` is the number of genes handed in; `n_tested` is how
#'   many of those could enter the test at all -- inputs present in the universe
#'   and in some set, which is the `DESize` that `enrichR`'s `GeneRatio` divides
#'   by; `n_sig_covered` is how many are covered by a *significant* regulator;
#'   `uncovered_sig` is `n_input - n_sig_covered`. `coverage_tested` is
#'   `n_tested / n_input`, not the larger ratio an "appears in any set" count
#'   would give.
#' @details `design` accepts `contrast`, `system`, `species`, `treatment`
#'   (`"none"` by default), `gene`, `direction`, `gene_col` and
#'   `coverage_floor`. `treatment` decides whether perturbation evidence (level
#'   A) is reachable; it is never inferred from free text. `background` names
#'   entities that belong to the setting rather than to the findings -- the
#'   receptor a drug acts through, a known positive control -- so that
#'   `unsupported_claims` does not report them as unevidenced. `coverage_floor` is
#'   `NULL` by default: a collection is only declared `untested` as a whole when
#'   the caller states a floor.
#' @seealso [mechanism_sets()] for the collections, [mechanism_xref()] for the
#'   level rules.
#' @examples
#' \dontrun{
#' ev <- evidence(c("HIF1A", "PGK1", "CA9", "BNIP3", "SLC2A1", "LDHA"))
#' ev$regulators
#' }
#' @export
evidence <- function(genes,
                     design       = list(),
                     pathway      = FALSE,
                     upstream     = "tf",
                     universe     = NULL,
                     network      = FALSE,
                     minGSSize    = 5,
                     maxGSSize    = Inf,
                     ## Which parts of the TF collection to use. `"gtrd"` keeps
                     ## only the ChIP-defined sets, whose identifiers are gene
                     ## symbols; the default uses the whole collection.
                     tft          = NULL,
                     print_pvalue = 1,
                     print_qvalue = 1) {
    upstream <- unique(upstream)
    if (!length(upstream)) {
        rlang::abort(c(
            "`upstream` must name at least one collection.",
            i = "Try `upstream = c(\"tf\", \"mirna\")`, the default."
        ))
    }
    ## "perturbagen" is measurement rather than a collection: it comes from
    ## GPSAdb over MCP and carries its own level rather than being graded by
    ## mechanism_xref().
    supported <- c("tf", "tfunion", "mirna", "perturbagen", "cellstate", "tissue")
    if (length(setdiff(upstream, supported))) {
        rlang::abort(c(
            paste0("Unknown `upstream`: ",
                   paste(setdiff(upstream, supported), collapse = ", "), "."),
            i = paste0("This release supports: ", paste(supported, collapse = ", "), ".")
        ))
    }

    ## `tft` selects part of C3:TFT and means nothing for the union, whose other
    ## four libraries have no motif/GTRD split. Silently ignoring it would answer
    ## a different question from the one asked: `tft = "gtrd"` is how a caller
    ## avoids the family and consensus-sequence labels that no perturbation
    ## dataset can ever answer for, and dropping that filter without a word
    ## turns a 504-set search into a 2869-set one.
    if (!is.null(tft) && "tfunion" %in% upstream) {
        rlang::abort(c(
            "`tft` cannot be combined with `upstream = \"tfunion\"`.",
            i = paste0("`tft` selects part of C3:TFT; the union also contains",
                       " ENCODE, ChEA, TF_Perturbations and TRRUST sets, which",
                       " have no GTRD/motif split."),
            i = paste0("Use `upstream = \"tf\"` with `tft`, or drop `tft` to",
                       " search the whole union.")
        ))
    }

    ## Whether the caller named the fold-change column themselves, before the
    ## defaults fill it in. An explicit name that is missing is an error; the
    ## default name being absent just means the table has no direction.
    logfc_given <- "logfc_col" %in% names(design)

    design <- mechanism_design(design)

    ## A named numeric vector is two things at once: the genes and the sign each
    ## one moved by. Held apart from `genes` because the rest of the function
    ## works in symbols, and the statistics are needed to split the list.
    signed <- NULL
    ## A data.frame carrying a fold-change column is the canonical input: it is
    ## the shape a DE table already has, and it holds the genes and the
    ## direction at once.
    if (is.data.frame(genes) && !is.null(design$logfc_col)) {
        has_col <- design$logfc_col %in% names(genes)
        if (!has_col && logfc_given) {
            rlang::abort(c(
                paste0("`genes` has no column \"", design$logfc_col, "\"."),
                i = "Set `design$logfc_col` to the fold-change column, or drop it to run without direction."
            ))
        }
        if (!has_col) {
            rlang::warn(paste0(
                "No \"", design$logfc_col, "\" column in `genes`, so the list has no ",
                "direction: the up/down split, the perturbation sign test and the ",
                "pathway layer cannot be asked. Set `design$logfc_col` to the ",
                "column that holds it."))
        } else {
            lfc <- suppressWarnings(as.numeric(genes[[design$logfc_col]]))
            sym <- as.character(genes[[design$gene_col]])
            keep <- !is.na(lfc) & !is.na(sym) & nzchar(sym)
            if (any(keep)) {
                signed <- stats::setNames(lfc[keep], sym[keep])
                signed <- signed[!duplicated(names(signed))]
            }
        }
    }
    if (is.numeric(genes) && !is.null(names(genes)) && any(nzchar(names(genes)))) {
        signed <- genes[!is.na(genes)]
        signed <- signed[nzchar(names(signed))]
        if (!length(signed)) signed <- NULL
    }

    ## An enrichResult carries its own gene list, its own background and its own
    ## keytype. All three are taken from the object rather than guessed, because
    ## getting any of them wrong changes the statistics silently: the default
    ## universe is every symbol in the annotation, which is not the background
    ## the enrichment itself was run against.
    from_result <- NULL
    ranked <- NULL
    if (inherits(genes, "gseaResult")) {
        ## A GSEA result already carries the question this function asks: every
        ## gene, ranked, with a statistic saying how far it moved. Taking the
        ## significant genes out of it would throw away exactly the information
        ## GSEA was run to keep.
        from_result <- .mechanism_from_gsearesult(genes, design)
        genes <- from_result$genes
        if (is.null(universe)) universe <- from_result$universe
        ## The ranking itself is what GSEA measured, so it is carried into the
        ## collection stage: re-running ORA on the core genes would throw away
        ## the direction and the continuous signal GSEA was run to obtain.
        ranked <- from_result$ranked
    } else if (inherits(genes, "enrichResult")) {
        ## Reading genes back out of an enrichment result loses them: @result is
        ## truncated at pvalueCutoff, so on the airway list 897 input genes came
        ## back as 483 with nothing in the output to say so. Taking the genes
        ## from the DE table and the pathways from the result loses nothing.
        rlang::abort(c(
            "`evidence()` no longer accepts an `enrichResult` as `genes`.",
            i = paste0("Pass the genes themselves and hand the result to `pathway =`, ",
                       "which uses it for the pathways while the genes stay complete: ",
                       "evidence(de, pathway = result)."),
            i = paste0("Reading genes back out of a result silently drops any gene that ",
                       "appears only in a pathway below its p-value cutoff.")
        ))
    }

    genes <- mechanism_genes(genes, design)
    universe <- mechanism_universe(universe, genes, design)

    ## The two halves of the caller's contrast, computed once. Both the
    ## perturbation sign test and the pathway layer ask about a half rather than
    ## about the list, and recomputing them separately invites them to differ.
    directed_up   <- if (!is.null(signed)) names(signed)[signed > 0] else NULL
    directed_down <- if (!is.null(signed)) names(signed)[signed < 0] else NULL

    needed <- setdiff(genes, universe)
    if (length(needed)) {
        rlang::warn(paste0(
            length(needed), " input gene(s) are outside `universe` and cannot be tested: ",
            paste(utils::head(needed, 5), collapse = ", "),
            if (length(needed) > 5) ", ..." else ""
        ))
    }

    ## Every source consulted below records its failures here, so an answer
    ## that could not be obtained is never reported as an answer that was empty.
    failures <- mechanism_failure_collector()

    collection_sources <- intersect(upstream,
                                    c("tf", "tfunion", "mirna", "cellstate", "tissue"))
    sets <- lapply(collection_sources, function(src) {
        mechanism_sets(src, species = design$species,
                       subset = if (identical(src, "tf")) tft else NULL)
    })
    names(sets) <- collection_sources

    mechanism_progress("evidence: testing collection(s): ",
                       paste(collection_sources, collapse = ", "),
                       " (", length(genes), " genes)")

    results <- lapply(collection_sources, function(src) {
        ## When the caller supplied a GSEA result, the collection is tested with
        ## GSEA too. Substituting ORA here was the first implementation and it
        ## was wrong: ORA takes the core genes as a binary list, so the ranking
        ## and the sign that GSEA produced are discarded, and the "GSEA path"
        ## becomes the ORA path with a different input. Measured both ways they
        ## agreed to within a rank, which is exactly what running the same test
        ## twice looks like.
        if (!is.null(ranked)) {
            return(.mechanism_gsea_collection(ranked, sets[[src]],
                                              minGSSize, maxGSSize))
        }
        ## Split by sign when the input carries one. Tested on the hypoxia
        ## contrast, this moved HIF1_Q3 from padj 1.2e-4 to 1.0e-5 -- a twelve-
        ## fold gain in power -- because a regulator's targets rise under
        ## hypoxia while other genes fall, and mixing the two directions asks a
        ## diluted question. It also supplies the direction ORA cannot otherwise
        ## give, without GSEA's ranking or permutations.
        if (!is.null(signed)) {
            return(.mechanism_directional_collection(signed, universe,
                                                     sets[[src]], minGSSize,
                                                     maxGSSize))
        }
        enricher(
            genes,
            universe     = universe,
            gson         = sets[[src]],
            minGSSize    = minGSSize,
            maxGSSize    = maxGSSize,
            ## Correction is applied to the whole collection; the print cutoffs
            ## are applied later, on the returned table, so that "not shown"
            ## never becomes "not tested".
            pvalueCutoff = 1,
            qvalueCutoff = 1
        )
    })
    names(results) <- collection_sources

    regulators <- do.call(rbind, lapply(collection_sources, function(src) {
        mechanism_regulator_table(results[[src]], sets[[src]], src, genes,
                                   species = design$species)
    }))
    ## No collection source leaves rbind() of nothing as NULL; the rest of the
    ## function treats `regulators` as a data.frame, so give it the empty one.
    if (is.null(regulators)) regulators <- mechanism_empty_regulators()
    ## biotrigger is available through evidence_biotrigger() but is not part of
    ## the bundle. Measured on the airway list, wiring it in raised `supported`
    ## from 21 to 25 and added no transcription factor -- TREM2, SORL1, PTPN23
    ## and STAG2 have no bearing on airway smooth muscle, and the two "sources"
    ## often cited the same dataset. A channel that raises a count without
    ## improving the mechanism is not an improvement.

    measured <- NULL
    if ("perturbagen" %in% upstream) {
        if (identical(design$treatment, "none")) {
            rlang::warn(paste0(
                "`upstream` includes \"perturbagen\" but `design$treatment` is \"none\". ",
                "Every perturbation dataset will be queried and graded A, although ",
                "the design does not state a treatment; set it, e.g. ",
                "design = list(treatment = \"knockdown\"), so the grade matches the study."
            ))
        }
        ## Asked about nominated regulators rather than a sample of the gene
        ## list: whether a candidate was ever knocked down is a well-posed
        ## question, while which genes happen to have datasets is a property of
        ## the database.
        cand <- mechanism_perturb_candidates(regulators)
        mechanism_progress("evidence: querying GPSAdb for ", length(cand),
                           " candidate regulator(s); this is the slow step")
        measured <- tryCatch(
            evidence_perturbation(
                genes,
                method     = mechanism_treatment_method(design),
                candidates = cand,
                considered = length(cand),
                system     = design$system %||% design$tissue,
                up         = directed_up,
                down       = directed_down
            ),
            error = function(e) {
                failures$add(source = "perturbagen", stage = "query",
                             message = conditionMessage(e))
                mechanism_empty_perturbation()
            }
        )
        if (nrow(failures$all())) {
            ## Say it loudly: the perturbation channel is the only one that can
            ## reach level A, so losing it changes what the run can conclude.
            rlang::warn(paste0(
                "Perturbation evidence could not be obtained: ",
                paste(unique(failures$all()$message), collapse = "; "),
                ". The bundle records this in $diagnostics$failures."))
        }
    }


    ## Nothing worked. An empty bundle here would report `absent` for every
    ## entity, which reads as a negative result when the truth is that no source
    ## answered.
    mechanism_abort_if_all_failed(
        failures, c(collection_sources,
                    intersect(upstream, c("perturbagen", "biotrigger"))))

    ## The canonical shape, whatever came out of rbind: `direction` belongs to
    ## the contract even when no measured row is present to populate it.
    regulators <- mechanism_pad_columns(regulators)
    if (!is.null(measured) && nrow(measured)) {
        regulators <- rbind(mechanism_pad_columns(regulators), measured[, names(mechanism_pad_columns(regulators)), drop = FALSE])
    }

    ## Network evidence is graded separately: no STRING channel measures a
    ## perturbation, so an edge can reach B and never A.
    ppi <- if (isTRUE(network)) {
        ev <- tryCatch(
            evidence_ppi(genes, species = mechanism_taxid(design$species)),
            error = function(e) e
        )
        if (inherits(ev, "error")) {
            failures$add(source = "ppi", stage = "network",
                         message = conditionMessage(ev))
            mechanism_empty_ppi()
        } else ev
    } else {
        NULL
    }

    out <- list(
        input = list(
            genes     = genes,
            universe  = universe,
            design    = design,
            namespace = "SYMBOL",
            ## Whether the caller's list carried a direction. The layers that
            ## need one report `untested` when it is FALSE, and print() has to be
            ## able to say which kind of silence this was.
            directed  = !is.null(signed)
        ),
        sets = lapply(sets, function(g) list(
            gsname   = g@gsname,
            version  = g@version,
            n_sets   = length(g@gsid2name$gsid),
            keytype  = g@keytype
        )),
        regulators  = regulators,
        xref        = mechanism_xref(),
        diagnostics = c(
            mechanism_diagnostics(regulators, names(sets), genes, sets, design, universe),
            ## Anything that could not be obtained, so a caller can tell an
            ## empty answer from an answer that was never received.
            list(failures = failures$all())
        ),
        provenance  = list(
            ## One entry per collection: they normally agree, but nothing
            ## stops a caller from mixing releases, and a scalar here would
            ## silently report only the first.
            ## No collection at all is legal -- a caller may ask for
            ## perturbation evidence alone -- so these must survive an empty
            ## `sets` list instead of calling vapply() on nothing.
            release       = if (length(sets)) vapply(sets, function(g) g@version, "") else character(0),
            release_consistent = length(unique(if (length(sets)) {
                vapply(sets, function(g) g@version, "")
            } else character(0))) <= 1L,
            accessed_date = unique(if (length(sets)) {
                vapply(sets, function(g) g@accessed_date, "")
            } else character(0)),
            upstream      = upstream,
            minGSSize     = minGSSize,
            maxGSSize     = maxGSSize,
            print_pvalue  = print_pvalue,
            print_qvalue  = print_qvalue,
            package       = "clusterProfiler",
            package_version = as.character(utils::packageVersion("clusterProfiler")),
            called_at     = as.character(Sys.time())
        )
    )
    for (src in collection_sources) out[[src]] <- results[[src]]
    if (!is.null(measured)) out$perturbagen <- measured
    ## Kept separate from the verdict table on purpose: this ranks experiments,
    ## not regulators, and a row here must never be read as a mechanism claim.
    if (!is.null(ppi)) out$ppi <- ppi

    ## The pathway layer, also separate. See evidence_pathway.R for why it must
    ## not become rows in $regulators.
    ## NULL means "not asked", the same as FALSE: a caller who computed a result
    ## conditionally ends up passing NULL, and treating that as a request would
    ## silently skip the layer while pretending to attempt it.
    if (!isFALSE(pathway) && !is.null(pathway)) {
        up <- directed_up %||% character(0)
        dn <- directed_down %||% character(0)
        if (isTRUE(pathway) && is.null(signed)) {
            rlang::warn(paste0(
                "`pathway = TRUE` needs a direction to enrich each pathway against; ",
                "the list has none, so the pathway layer was skipped. Supply ",
                "`design$logfc_col`, or pass a computed result as `pathway =`."))
        } else {
            pw <- tryCatch({
                pathways <- .mechanism_pathways(pathway, up, dn, universe, design)
                .mechanism_pathway_layer(pathways, universe, sets,
                                         minGSSize, maxGSSize,
                                         top = design$pathway_top)
            }, error = function(e) {
                failures$add(source = "pathway", stage = "enrich",
                             message = conditionMessage(e))
                NULL
            })
            if (!is.null(pw)) out$pathways <- pw
        }
    }

    ## The class must exist before the writable table is built: it defers to
    ## consensus(), which refuses anything that is not an evidenceBundle.
    class(out) <- c("evidenceBundle", "list")
    out$writable <- .evidence_writable_table(out)
    out
}

#' Print an evidence bundle
#'
#' Show the input size, the collections used and the most significant
#' regulators, one row per regulator with its evidence level. The cutoffs come
#' from `print_pvalue` / `print_qvalue` recorded when the bundle was built;
#' they affect this display only, never the table stored in `x$regulators`.
#'
#' @param x An `evidenceBundle` from [evidence()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.evidenceBundle <- function(x, ...) {
    d <- x$diagnostics$counts
    cat("<evidenceBundle>\n")
    ## n_input is the same for every source row; [1] collapses that repetition.
    ## With no collection source -- perturbation evidence only -- there is no
    ## row, so fall back to the input itself.
    n_in <- if (!is.null(d) && nrow(d)) d$n_input[1] else length(x$input$genes)
    cat("  input genes:  ", n_in, "  universe: ", length(x$input$universe), "\n")
    rel <- unique(x$provenance$release)
    rel <- if (length(rel)) paste0(" (", paste(rel, collapse = ", "), ")") else ""
    cat("  collections:  ", paste(x$provenance$upstream, collapse = ", "), rel, "\n", sep = "")

    r <- x$regulators
    noted <- sum(!is.na(r$padj_within_source))
    shown <- r[!is.na(r$padj_within_source) &
                   r$padj_within_source <= (x$provenance$print_pvalue %||% 1), , drop = FALSE]
    ## "N of N" said nothing: only sets that were actually tested carry a
    ## p-value, so the tested count is the informative denominator.
    cat("  regulators:   ", nrow(shown), " shown of ", noted, " tested, ",
        sum(shown$named), " naming a factor\n", sep = "")

    ## A source that failed is not a source that found nothing.
    fails <- x$diagnostics$failures
    if (!is.null(fails) && nrow(fails)) {
        cat("  FAILED:       ", nrow(fails), " source(s) could not be consulted\n", sep = "")
        for (i in seq_len(nrow(fails))) {
            cat("     ", fails$source[i], " [", fails$stage[i], "]: ",
                substr(fails$message[i], 1, 60), "\n", sep = "")
        }
    }

    ## A similarity ranking answers a different question from everything above,
    ## so it is reported apart from the verdicts.
    ## How much of the measured channel had data to answer with. Three hits out
    ## of three hundred candidates is not the same claim as three out of three,
    ## and nothing else in this output would tell the two apart.
    if (!is.null(x$perturbagen)) {
        q <- attr(x$perturbagen, "queried")
        cons <- attr(x$perturbagen, "considered")
        if (!is.null(q) && !is.na(q)) {
            cat("  measured:     ", nrow(x$perturbagen),
                " dataset(s) for ", q, " of ", cons,
                " candidate regulator(s) queried\n", sep = "")
        }
    }

    ## What the optional layers could and could not be asked. A layer nobody
    ## asked is not a layer that found nothing, and without these two lines the
    ## reader cannot tell the difference.
    if (!is.null(x$pathways)) {
        cat("  pathways:     ", length(unique(x$pathways$pathway)),
            " pathway(s) broken down into candidate regulators\n", sep = "")
    }
    if (!isTRUE(x$input$directed)) {
        cat("  NO DIRECTION:  the up/down split, the perturbation sign test and",
            "\n                 the pathway layer could not be asked.\n", sep = "")
    }

    if (nrow(shown)) {
        shown <- shown[order(shown$padj_within_source), , drop = FALSE]
        shown <- utils::head(shown, 10)
        print(
            data.frame(
                regulator = shown$regulator,
                source    = shown$source,
                level     = shown$level,
                named     = shown$named,
                padj      = signif(shown$padj_within_source, 3),
                Count     = shown$Count,
                setSize   = shown$setSize,
                stringsAsFactors = FALSE
            ),
            row.names = FALSE
        )
    }
    invisible(x)
}

## design: fill the defaults once, so every later step sees the same shape.
mechanism_design <- function(design) {
    defaults <- list(
        contrast       = NULL,
        system         = NULL,
        species        = "Homo sapiens",
        treatment      = "none",
        gene           = NULL,
        direction      = NULL,
        gene_col       = "gene",
        ## Which column of a data.frame carries the direction. Named rather than
        ## guessed: a table that spells it "logFC" would otherwise run with every
        ## directional layer reporting `untested`, which reads as "not
        ## answerable" when the truth is "a column was not found".
        logfc_col      = "log2FoldChange",
        ## The conventional enrichment this function runs when `pathway = TRUE`.
        pathway_ontology = "BP",
        pathway_top    = 20L,
        coverage_floor = NULL,
        ## Entities that belong to the experimental setting rather than to the
        ## findings, e.g. the receptor a treatment acts through.
        background     = character(0)
    )
    unknown <- setdiff(names(design), names(defaults))
    if (length(unknown)) {
        rlang::warn(paste0("Ignoring unknown `design` field(s): ",
                           paste(unknown, collapse = ", ")))
    }
    utils::modifyList(defaults, design[intersect(names(design), names(defaults))])
}

mechanism_genes <- function(genes, design) {
    genes <- if (is.numeric(genes) && !is.null(names(genes))) {
        ## The names are the genes; the values are the statistics, which the
        ## caller's `universe` and the directional split use.
        names(genes)
    } else if (inherits(genes, "gseaResult")) {
        ## Only reached when the caller passed a result object but asked for no
        ## ranking: the ranked list is the input GSEA itself used.
        names(genes@geneList)
    } else if (inherits(genes, "enrichResult")) {
        unique(unlist(strsplit(as.character(genes$geneID), "/")))
    } else if (is.data.frame(genes)) {
        col <- design$gene_col
        if (!col %in% names(genes)) {
            rlang::abort(c(
                paste0("`genes` has no column \"", col, "\"."),
                i = "Set `design$gene_col` to the column holding the gene symbols."
            ))
        }
        as.character(genes[[col]])
    } else {
        as.character(genes)
    }
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    if (!length(genes)) rlang::abort("`genes` is empty after removing NA and blank entries.")
    genes
}

## The universe is returned as given (or as every symbol in org.Hs.eg.db).
##
## It is deliberately NOT unioned with `genes`: adding an input gene that the
## caller excluded from the background would enlarge the denominator of the
## hypergeometric test, which is the one thing a universe exists to pin down.
## Input genes outside the universe are dropped by enricher() itself, and
## evidence() warns about them.
mechanism_universe <- function(universe, genes, design) {
    if (is.null(universe)) {
        pkg <- "org.Hs.eg.db"
        if (!requireNamespace(pkg, quietly = TRUE)) {
            rlang::abort(c(
                "No `universe` supplied and org.Hs.eg.db is not installed.",
                i = "Install it with BiocManager::install(\"org.Hs.eg.db\"), or pass `universe`."
            ))
        }
        universe <- AnnotationDbi::keys(getExportedValue(pkg, pkg), keytype = "SYMBOL")
    }
    universe <- as.character(universe)
    universe <- unique(universe[!is.na(universe) & nzchar(universe)])
    if (!length(universe)) rlang::abort("`universe` is empty after removing NA and blank entries.")
    universe
}

## One row per tested set, with the regulator label it collapses to.
mechanism_regulator_table <- function(res, gson_obj, source, selected,
                                     species = "Homo sapiens") {
    ## A collection in which no set passes the size window makes enricher()
    ## return NULL; that is an empty answer, not an error, so it becomes an
    ## empty table and the diagnostics still report what was tested.
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
        return(data.frame(
            regulator = character(), source = character(), term = character(),
            named = logical(), level = character(),
            padj_within_source = numeric(), setSize = integer(), Count = integer(),
            family = logical(), direction = character(), stringsAsFactors = FALSE
        ))
    }
    df <- as.data.frame(res)
    df <- df[!is.na(df$ID), , drop = FALSE]
    members <- gson_obj@gsid2gene

    ## A gseaResult has no `Count`: the number of input genes in the set has to
    ## come from the leading edge, and the direction from the sign of NES.
    is_gsea <- inherits(res, "gseaResult")
    if (is_gsea) {
        df$Count <- vapply(strsplit(as.character(df$core_enrichment), "/"),
                           length, integer(1))
        df$direction <- .mechanism_gsea_direction(df$NES)
    } else if (!"direction" %in% names(df)) {
        ## The directional path already knows which half of the list a set was
        ## enriched in; an undirected ORA has no sign to report.
        df$direction <- NA_character_
    }

    ## setSize counts unique (gsid, gene) pairs; gson::gson() already dedupes
    ## gsid2gene, and this keeps the count correct if that ever stops.
    set_size <- table(unique(members[, c("gsid", "gene")])$gsid)
    ## TF sets need collapsing; cell-state and tissue identifiers already are
    ## the entity name, and only a TF label can be "unnamed" in the sense of
    ## naming no factor at all.
    regulator <- mechanism_entity(df$ID, source)
    unnamed <- if (identical(source, "tf")) mechanism_unnamed(regulator) else rep(FALSE, length(regulator))

    data.frame(
        regulator          = mechanism_alias_regulator(regulator),
        source             = source,
        term               = df$ID,
        named              = !unnamed,
        ## A family or complex label cannot be looked up in a perturbation
        ## database, and its string must not be matched against gene symbols.
        family             = mechanism_is_family_label(regulator, species),
        ## Whether this row is a measurement, a prediction or a signature. The
        ## level says how strong it is and the source says which collection
        ## produced it; neither says whether anything was measured.
        evidence_kind      = mechanism_kind_of(source, df$ID),
        level              = mechanism_reported_level(source, df$ID, regulator),
        padj_within_source = df$p.adjust,
        setSize            = as.integer(set_size[df$ID]),
        Count              = df$Count,
        ## Filled only on the GSEA path; ORA has no sign to report.
        direction          = df$direction,
        stringsAsFactors   = FALSE
    )
}

## The counts whose denominators differ. Kept in one data.frame so a caller
## cannot pick a numerator and a denominator from different questions.
mechanism_diagnostics <- function(regulators, upstream, genes, sets, design, universe) {
    n_input <- length(genes)

    ## Input genes the universe keeps. A gene outside the universe never
    ## entered the test, whatever the collection says about it, and this is
    ## the same intersection enricher() applies before computing DESize.
    testable <- intersect(genes, universe)

    per_source <- lapply(upstream, function(src) {
        members <- unique(sets[[src]]@gsid2gene$gene)
        if (is.null(members)) members <- character(0)
        r <- regulators[regulators$source == src, , drop = FALSE]
        sig <- r[mechanism_significant(r), , drop = FALSE]

        covered_by_set <- intersect(testable, members)
        covered_sig <- if (nrow(sig)) {
            intersect(testable, unique(sets[[src]]@gsid2gene$gene[
                sets[[src]]@gsid2gene$gsid %in% sig$term]))
        } else {
            character(0)
        }
        ## DESize: the number of input genes that could enter the test at all,
        ## i.e. inputs found in the universe *and* in some set. This is the
        ## denominator enrichit puts in GeneRatio, and it is the only honest
        ## measure of "how much of the list was tested" -- counting genes the
        ## universe excluded reports a coverage the test never had.
        n_tested <- length(covered_by_set)
        annotated <- n_tested / n_input
        floor <- design$coverage_floor
        status <- if (is.null(floor)) {
            "graded"
        } else if (annotated >= floor) {
            "graded"
        } else {
            "untested"
        }

        data.frame(
            source           = src,
            n_input          = n_input,
            n_tested         = n_tested,
            n_sig_covered    = length(covered_sig),
            uncovered_sig    = n_input - length(covered_sig),
            coverage_tested  = round(annotated, 4),
            coverage_floor   = if (is.null(floor)) NA_real_ else floor,
            status           = status,
            stringsAsFactors = FALSE
        )
    })

    list(
        counts   = mechanism_counts_or_empty(per_source),
        ## Where each individual gene stands, for the chapter that has to
        ## explain "absent is evidence, untested is not".
        per_gene = mechanism_per_gene(regulators, upstream, genes, sets, design, universe)
    )
}

mechanism_per_gene <- function(regulators, upstream, genes, sets, design, universe) {
    rows <- lapply(upstream, function(src) {
        members_of <- sets[[src]]@gsid2gene
        r <- regulators[regulators$source == src, , drop = FALSE]
        sig <- r[mechanism_significant(r), , drop = FALSE]
        sig_genes <- unique(members_of$gene[members_of$gsid %in% sig$term])
        ## A gene outside the universe is untested even if a set lists it.
        annotated <- genes %in% universe & genes %in% unique(members_of$gene)

        data.frame(
            gene   = genes,
            source = src,
            status = ifelse(
                genes %in% sig_genes, "supported",
                ifelse(annotated, "absent", "untested")
            ),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

## STRING wants a taxonomy id; design$species carries the organism name.
mechanism_taxid <- function(species) {
    if (is.null(species) || identical(species, "auto")) return("auto")
    if (grepl("^[0-9]+$", species)) return(species)
    if (requireNamespace("yulab.utils", quietly = TRUE)) {
        id <- tryCatch(getTaxID(species), error = function(e) NA)
        if (!is.na(id)) return(as.character(id))
    }
    "auto"
}

## The sentence table is consensus()'s, not a second opinion computed here.
##
## Grading each term row separately produced contradictory advice for one
## regulator -- 'AP4' was simultaneously "state it as a hypothesis" (one
## significant set) and "no supporting evidence" (three non-significant ones),
## because a verdict computed per row cannot see the entity's other rows.
## consensus() decides per entity and carries its own `writable` column, so the
## bundle simply defers to it.
.evidence_writable_table <- function(bundle) {
    cv <- consensus(bundle)
    if (!nrow(cv)) {
        return(data.frame(regulator = character(), sources = character(),
                          n_sources = integer(), max_level = character(),
                          verdict = character(), padj = numeric(),
                          direction = character(), writable = character(),
                          stringsAsFactors = FALSE))
    }
    cv
}

## A frame with every column the regulators table uses, so measured rows and
## collection rows can be bound together without a shared subset.
mechanism_pad_columns <- function(r) {
    template <- data.frame(
        regulator = character(), source = character(), term = character(),
        named = logical(), level = character(), padj_within_source = numeric(),
        setSize = integer(), Count = integer(),
        ## A family or complex label can never be looked up as a gene, so this
        ## column travels with every row and not only the collection ones.
        family = logical(),
        evidence_kind = character(),
        ## Perturbation rows carry the direction the responding genes moved;
        ## collection rows have none. The column must survive the rbind, so it
        ## belongs in the template and not only in the measured table.
        direction = character(),
        ## Where a measured row came from: a dataset id means nothing to a
        ## reader without the cell line and the method behind it.
        dataset = character(), system = character(), method = character(),
        ## Sources that report a score rather than a test say for themselves
        ## whether a row is a hit; see mechanism_significant().
        significant = logical(), effect_size = numeric(),
        ## Whether a measured row was run in the system the caller named.
        ## NA when no system was named, so the column never claims a match
        ## that was not checked.
        in_system = logical(),
        ## The two halves' sign tests; collection rows carry NA.
        up_rise = integer(), up_fall = integer(), up_sign_padj = numeric(),
        down_rise = integer(), down_fall = integer(), down_sign_padj = numeric(),
        stringsAsFactors = FALSE
    )
    ## Assign the right length, not the template's zero: `r` may already have
    ## rows, and a zero-length vector cannot be recycled into them.
    for (col in setdiff(names(template), names(r))) {
        r[[col]] <- rep(template[[col]], length.out = nrow(r))
    }
    r[, names(template), drop = FALSE]
}

## One row per source is the contract. With no collection there is no row, but a
## NULL breaks callers that subset counts$source, so return the empty frame.
mechanism_counts_or_empty <- function(per_source) {
    counts <- do.call(rbind, per_source)
    if (!is.null(counts)) return(counts)
    data.frame(
        source = character(), n_input = integer(), n_tested = integer(),
        n_sig_covered = integer(), uncovered_sig = integer(),
        coverage_tested = numeric(), coverage_floor = numeric(),
        status = character(), stringsAsFactors = FALSE
    )
}

## design$treatment names the perturbation kind; GPSAdb wants a method string.
## "none" means the study states no treatment, so nothing is queried.
mechanism_treatment_method <- function(design) {
    tr <- design$treatment %||% "none"
    if (identical(tr, "none")) return(NULL)
    tr
}

## Regulators worth asking GPSAdb about.
##
## Only a regulator the collections already found *significant* is worth a
## query: asking about a factor whose target sets were unremarkable adds a
## request and cannot strengthen the result. Candidates are ordered by
## significance so the cap keeps the strongest, and the cap exists because each
## candidate costs one request -- asking about every named regulator in a
## 900-gene list meant 613 requests and a run that never finished.
##
## shortcut: a capped list, not the full set. Raise
## `options(clusterProfiler_perturb_candidates =)` when completeness matters
## more than the wait; how many were queried and how many were skipped is
## recorded on the returned table.
mechanism_perturb_candidates <- function(regulators, alpha = 0.05,
    max_n = getOption("clusterProfiler_perturb_candidates", 30L)) {
    if (is.null(regulators) || !nrow(regulators)) return(NULL)
    ## A family or complex label is never worth a query: GPSAdb indexes gene
    ## symbols, so a request for 'AP1' or 'NFAT' cannot return anything.
    keep <- regulators$named &
        !(regulators$family %||% FALSE) &
        !is.na(regulators$padj_within_source) &
        regulators$padj_within_source < alpha
    if (!any(keep)) return(NULL)
    r <- regulators[keep, , drop = FALSE]
    r <- r[order(r$padj_within_source), , drop = FALSE]
    named <- unique(r$regulator[!is.na(r$regulator) & nzchar(r$regulator)])
    if (!length(named)) return(NULL)
    considered <- length(named)
    if (considered > max_n) {
        rlang::warn(paste0(
            considered, " regulators have significant target sets; querying GPSAdb ",
            "for the ", max_n, " strongest. Set ",
            "options(clusterProfiler_perturb_candidates = ", considered, ") ",
            "to query them all."))
        named <- named[seq_len(max_n)]
    }
    ## How many were considered is what makes the result interpretable: what the
    ## sweep covered is part of the answer, not an implementation detail.
    attr(named, "considered") <- considered
    named
}

## Genes, background and keytype taken from an enrichment result.
##
## The obvious source is wrong. `@gene` is `unique(c(original_genes, every
## member of every gene set, the tested genes))`, so on a ten-gene input it
## holds 15941 entries -- it is the annotation, not the query. The recoverable
## source is the `geneID` column of `@result`, which lists, per pathway, the
## input genes that pathway contains.
##
## That source is lossy, and the loss is silent. `@result` is truncated by
## `pvalueCutoff` at construction, so a gene appearing only in pathways that did
## not reach the cutoff cannot be recovered. Measured on the airway list with
## the default cutoff of 0.05, 897 input genes came back as 483 -- 414 lost,
## with nothing in the output to say so. Rebuilding the result with
## `pvalueCutoff = 1` recovered 723. The warning below exists because that
## difference changes the answer and the caller cannot see it.
.mechanism_from_enrichresult <- function(x, design) {
    df <- as.data.frame(x)
    if (!"geneID" %in% names(df)) {
        rlang::abort(c(
            "This enrichment result has no `geneID` column, so its genes cannot be read.",
            i = "Pass the gene list directly instead."
        ))
    }
    genes <- unique(unlist(strsplit(as.character(df$geneID), "/"), use.names = FALSE))
    genes <- genes[!is.na(genes) & nzchar(genes)]

    cutoff <- tryCatch(x@pvalueCutoff, error = function(e) 1)
    if (is.numeric(cutoff) && length(cutoff) == 1 && cutoff < 1) {
        rlang::warn(c(
            paste0("The enrichment result was built with pvalueCutoff = ", cutoff,
                   ", so only pathways that reached it are stored."),
            i = paste0("Only the ", length(genes),
                       " gene(s) appearing in those pathways can be recovered; ",
                       "any gene unique to a non-significant pathway is missing."),
            i = "Rebuild with pvalueCutoff = 1, qvalueCutoff = 1 for the full list."
        ))
    }

    keytype <- tryCatch(x@keytype, error = function(e) NULL)
    species <- tryCatch(x@organism, error = function(e) NULL)
    if (is.null(species) || !length(species) || is.na(species)) {
        species <- design$species %||% "Homo sapiens"
    }
    universe <- tryCatch(as.character(x@universe), error = function(e) character())

    ## evidence() works in symbols and its collections are symbol-keyed, so an
    ## ENTREZID result would otherwise be matched against nothing and report a
    ## confident `absent` for every entity.
    if (!is.null(keytype) && length(keytype) && nzchar(keytype) &&
        !identical(keytype, "SYMBOL")) {
        orgdb <- mechanism_orgdb(species)
        if (is.null(orgdb)) {
            rlang::abort(c(
                paste0("The result is keyed by ", keytype,
                       " and the annotation for ", species, " is not installed."),
                i = "Install the matching org.*.eg.db package, or pass symbols."
            ))
        }
        genes <- .mechanism_to_symbol(orgdb, genes, keytype)
        if (length(universe)) universe <- .mechanism_to_symbol(orgdb, universe, keytype)
    }

    list(genes = genes, universe = universe, keytype = keytype, species = species)
}

## Map one keytype to symbols, dropping what has no symbol rather than keeping
## an identifier that would be silently matched against nothing.
.mechanism_to_symbol <- function(orgdb, keys, keytype) {
    mapped <- tryCatch(
        AnnotationDbi::mapIds(orgdb, keys = keys, column = "SYMBOL",
                              keytype = keytype, multiVals = "first"),
        error = function(e) NULL)
    if (is.null(mapped)) {
        rlang::abort(paste0("Could not map `", keytype, "` to SYMBOL."))
    }
    mapped <- unique(as.character(mapped[!is.na(mapped)]))
    mapped[nzchar(mapped)]
}

## The ranked list, background and keytype taken from a GSEA result.
##
## Unlike enrichResult, this object was built from the whole ranked list, and it
## still holds it: `@geneList` is every gene with its statistic, sorted, which
## is precisely the input the mechanism analysis wants. The significant genes
## are not extracted, because extracting them would re-impose the cutoff that
## GSEA exists to avoid -- a gene that moved smoothly with the list belongs in
## the analysis even when it did not clear a p-value of its own.
##
## `@geneSets` holds every set tested, so the background is the ranked list's
## own gene universe rather than the entire annotation.
.mechanism_from_gsearesult <- function(x, design) {
    ranked <- tryCatch(x@geneList, error = function(e) NULL)
    if (is.null(ranked) || !length(ranked)) {
        rlang::abort(c(
            "This GSEA result carries no ranked list.",
            i = "Pass the named gene list directly instead."
        ))
    }
    genes <- names(ranked)
    if (is.null(genes) || all(is.na(genes)) || !any(nzchar(genes))) {
        rlang::abort(c(
            "The GSEA result's ranked list has no gene names.",
            i = "Pass a named list, e.g. `setNames(de$stat, de$gene)`."
        ))
    }
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])

    keytype <- tryCatch(x@keytype, error = function(e) NULL)
    species <- tryCatch(x@organism, error = function(e) NULL)
    if (is.null(species) || !length(species) || is.na(species)) {
        species <- design$species %||% "Homo sapiens"
    }

    ## The universe is the ranked list itself. GSEA tested every set against
    ## these genes, so this is the background its own p-values were computed
    ## against; substituting the annotation's full symbol list would change the
    ## denominator the same way it did for enrichResult.
    universe <- genes

    ## But which genes are *input* to the mechanism question is a separate
    ## decision. The ranked list holds every gene, including the ten thousand
    ## that barely moved, and asking "which regulator explains all 22899" has
    ## no answer -- measured on the airway list it returned 0 supported
    ## regulators, because no factor's target set stands out against the whole
    ## transcriptome.
    ##
    ## `core_enrichment` is the alternative the object also provides: the genes
    ## that produced each enriched set's ranking, on 4594 of 22899 for this
    ## list. Those are the genes that moved together with a coherent programme,
    ## which is what a regulator is supposed to explain.
    core <- .gsea_core_genes(x)
    if (length(core)) genes <- core

    if (!is.null(keytype) && length(keytype) && nzchar(keytype) &&
        !identical(keytype, "SYMBOL")) {
        orgdb <- mechanism_orgdb(species)
        if (is.null(orgdb)) {
            rlang::abort(c(
                paste0("The result is keyed by ", keytype,
                       " and the annotation for ", species, " is not installed."),
                i = "Install the matching org.*.eg.db package, or pass symbols."
            ))
        }
        ## The statistic travels with the gene through the mapping, otherwise
        ## every renamed gene would silently drop out of the ranking.
        mapped <- .mechanism_to_symbol(orgdb, names(ranked), keytype)
        stats <- ranked[match(mapped, names(ranked))]
        ranked <- sort(stats[!is.na(stats)], decreasing = TRUE)
        genes <- unique(names(ranked))
        universe <- genes
    }

    list(genes = genes, universe = universe, keytype = keytype,
         species = species,
         ## The ranking, which is what distinguishes the GSEA path from the ORA
         ## path. Dropped here it silently degraded the whole function to ORA on
         ## the core genes, and the two "different" front ends then agreed to
         ## within a rank -- which is what running the same test twice looks
         ## like.
         ranked = ranked)
}

## Genes that produced the leading edge of the enriched sets.
##
## A GSEA result reports, per set, the genes that made it rank as it did. Their
## union across the significant sets is the moving part of the transcriptome --
## the analogue of the significant-gene list an ORA analysis starts from, but
## chosen by trend rather than by a threshold.
.gsea_core_genes <- function(x) {
    df <- tryCatch(as.data.frame(x), error = function(e) NULL)
    if (is.null(df) || !"core_enrichment" %in% names(df)) return(character(0))
    ## Only what the caller kept: a GSEA result truncated at pvalueCutoff is
    ## the caller's own statement of which sets matter.
    genes <- unique(unlist(strsplit(as.character(df$core_enrichment), "/"),
                          use.names = FALSE))
    genes[!is.na(genes) & nzchar(genes)]
}

## Test one collection with GSEA rather than ORA.
##
## Returns the same `enricher()`-shaped object the rest of the function expects,
## so the parsing below is untouched. `pvalueCutoff`/`qvalueCutoff` are set to 1
## for the same reason as the ORA call: "not shown" must never become "not
## tested".
.mechanism_gsea_collection <- function(ranked, gson, minGSSize, maxGSSize) {
    tryCatch(
        suppressWarnings(
            GSEA(ranked, gson = gson, minGSSize = minGSSize,
                 maxGSSize = maxGSSize, pvalueCutoff = 1, verbose = FALSE)
        ),
        error = function(e) NULL
    )
}

## A GSEA result's per-set direction.
##
## ORA cannot supply one: a binary gene list has no sign. GSEA's normalised
## enrichment score does -- a set whose members sit at the top of the ranking is
## up, at the bottom is down.
.mechanism_gsea_direction <- function(nes) {
    ifelse(is.na(nes), NA_character_,
           ifelse(nes > 0, "up", ifelse(nes < 0, "down", NA_character_)))
}

## Test a collection against an up list and a down list, separately.
##
## Three decisions, each of which changes the answer if made differently:
##
## * Both halves are corrected together. Splitting the list doubles the number
##   of hypotheses tested, so correcting each half on its own would report an
##   FDR for a question asked twice as often as the count suggests.
## * A set enriched in both directions keeps the stronger half. One contrast
##   cannot both activate and repress the same programme, and leaving both rows
##   makes consensus() report `mixed` for want of a rule rather than because the
##   biology is ambiguous.
## * Genes with a statistic of exactly zero are in neither half. They are rare,
##   and assigning them to one side would invent a direction.
.mechanism_directional_collection <- function(signed, universe, gson,
                                              minGSSize, maxGSSize) {
    halves <- list(list("up", names(signed)[signed > 0]),
                   list("down", names(signed)[signed < 0]))
    parts <- list()
    for (h in halves) {
        if (length(h[[2]]) < 2L) next
        r <- tryCatch(
            suppressWarnings(enricher(h[[2]], universe = universe, gson = gson,
                                      minGSSize = minGSSize, maxGSSize = maxGSSize,
                                      pvalueCutoff = 1, qvalueCutoff = 1)),
            error = function(e) NULL)
        if (is.null(r) || !nrow(as.data.frame(r))) next
        df <- as.data.frame(r)
        df$direction <- h[[1]]
        parts[[length(parts) + 1L]] <- df
    }
    if (!length(parts)) return(NULL)

    df <- do.call(rbind, parts)
    rownames(df) <- NULL
    df$p.adjust <- stats::p.adjust(df$pvalue, method = "BH")
    df <- df[order(df$p.adjust), , drop = FALSE]
    df[!duplicated(df$ID), , drop = FALSE]
}
