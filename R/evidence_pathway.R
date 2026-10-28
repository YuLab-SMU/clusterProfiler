## Pathway level: which regulator explains each enriched pathway.
##
## A different question from the collection layer. There, an upstream set is
## tested against the caller's whole gene list -- "whose target set is
## over-represented here". Here it is tested against the genes of one pathway --
## "does this regulator's target set explain this pathway". The statistic is the
## same, so this is enricher() with a different input and no second engine.
##
## The rows live in their own table and are never appended to $regulators. A
## pathway set is derived from the very collections the collection layer was
## built from, so counting it as a source would let one collection corroborate
## itself -- the failure that made evidence_biotrigger() raise the `supported`
## count without adding a single regulator.

## Which pathways to break down, and which way each one moved.
##
## `pathway = TRUE` runs the enrichment here, on the up and the down half
## separately: that is the directional form of the ordinary enrichment a caller
## would have run anyway, and running it on both halves is what lets a pathway
## carry a direction at all. An existing result is reused unchanged, so a caller
## who has already chosen their thresholds does not have them overridden.
.mechanism_pathways <- function(pathway, up, down, universe, design) {
    df <- if (isTRUE(pathway)) {
        .mechanism_pathways_compute(up, down, universe, design)
    } else if (inherits(pathway, c("enrichResult", "gseaResult"))) {
        .mechanism_pathways_from_result(pathway)
    } else {
        NULL
    }
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- .mechanism_pathways_direction(df, up, down)
    .mechanism_pathways_dedupe(df)
}

.mechanism_pathways_compute <- function(up, down, universe, design) {
    orgdb <- mechanism_orgdb(design$species)
    if (is.null(orgdb)) {
        rlang::abort(c(
            paste0("`pathway = TRUE` needs the annotation for ",
                   design$species, "."),
            i = "Install the matching org.*.eg.db package, or pass a computed result as `pathway =`."
        ))
    }
    ont <- design$pathway_ontology %||% "BP"
    parts <- list()
    for (h in list(list("up", up), list("down", down))) {
        genes <- h[[2]]
        ## enrichGO on a handful of genes returns nothing usable and costs a
        ## full annotation load; a half this small is reported as absent rather
        ## than queried.
        if (length(genes) < 5L) next
        r <- tryCatch(
            suppressWarnings(enrichGO(
                gene          = genes,
                OrgDb         = orgdb,
                keyType       = "SYMBOL",
                ont           = ont,
                universe      = universe,
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = FALSE)),
            error = function(e) NULL)
        if (is.null(r)) next
        d <- as.data.frame(r)
        if (!nrow(d)) next
        d$direction <- h[[1]]
        parts[[length(parts) + 1L]] <- d
    }
    if (!length(parts)) return(NULL)
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}

## Fill in a direction the source object does not carry.
##
## An `enrichResult` says nothing about which way its pathways moved -- there is
## no NES to read -- so a reused result would reach the narrative as "direction
## not stated" even when the caller handed in a signed gene list. The halves are
## already in hand, so the pathway's own genes say it: whichever half they fall
## in more often is the direction.
##
## A tie is left NA rather than broken. A pathway split evenly between the up
## and down halves moved neither way, and picking the larger count would invent
## one.
.mechanism_pathways_direction <- function(df, up, down) {
    if ("direction" %in% names(df) && any(!is.na(df$direction))) return(df)
    if (!length(up) && !length(down)) {
        df$direction <- NA_character_
        return(df)
    }
    df$direction <- vapply(strsplit(as.character(df$geneID), "/"), function(genes) {
        n_up <- sum(genes %in% up)
        n_down <- sum(genes %in% down)
        if (n_up > n_down) "up" else if (n_down > n_up) "down" else NA_character_
    }, "")
    df
}

## A result the caller already built. The direction comes from the sign of NES
## when the object has one, because for a GSEA result that is the only statement
## about direction the object carries.
.mechanism_pathways_from_result <- function(x) {
    d <- tryCatch(as.data.frame(x), error = function(e) NULL)
    if (is.null(d) || !nrow(d) || !all(c("ID", "geneID") %in% names(d))) {
        rlang::abort(c(
            "`pathway =` needs an enrichment result with `ID` and `geneID` columns.",
            i = "Pass the gene list with `pathway = TRUE` to compute the enrichment instead."
        ))
    }
    d$direction <- if ("NES" %in% names(d)) {
        .mechanism_gsea_direction(d$NES)
    } else {
        NA_character_
    }

    ## `enrichGO(readable = FALSE)` writes Entrez ids into `geneID`, and the
    ## collections are keyed by symbol, so passing those ids straight through
    ## matches nothing and the pathway layer comes back empty -- which reads as
    ## "no regulator explains any pathway" rather than "the ids were not
    ## translated".
    keytype <- tryCatch(x@keytype, error = function(e) NULL)
    if (!is.null(keytype) && length(keytype) && nzchar(keytype) &&
        !identical(keytype, "SYMBOL")) {
        species <- tryCatch(x@organism, error = function(e) NULL)
        orgdb <- mechanism_orgdb(if (length(species) && !is.na(species)) species
                                 else "Homo sapiens")
        if (is.null(orgdb)) {
            rlang::abort(c(
                paste0("`pathway =` was given a result keyed by ", keytype,
                       " and the annotation to translate it is not installed."),
                i = "Install the matching org.*.eg.db package, or pass `pathway = TRUE`."
            ))
        }
        d$geneID <- vapply(strsplit(as.character(d$geneID), "/"), function(ids) {
            sym <- tryCatch(
                AnnotationDbi::mapIds(orgdb, keys = ids, column = "SYMBOL",
                                      keytype = keytype, multiVals = "first"),
                error = function(e) NULL)
            sym <- if (is.null(sym)) character(0) else
                unique(as.character(sym[!is.na(sym)]))
            paste(sym[nzchar(sym)], collapse = "/")
        }, "")
        d <- d[nzchar(d$geneID), , drop = FALSE]
    }
    d
}

## Collapse pathways that hold exactly the same genes.
##
## GO's near-synonyms are the reason: on the hypoxia list 125 terms reduced to
## 84 distinct gene sets, and "response to hypoxia", "response to decreased
## oxygen levels" and "response to oxygen levels" are three names for one set.
## Testing each name would run the same hypergeometric three times and then
## report three findings that are one finding. The best-corrected name stands
## for the set, and the others are kept in `aliases` rather than dropped -- a
## reader searching for a term needs to find it.
.mechanism_pathways_dedupe <- function(df) {
    key <- vapply(strsplit(as.character(df$geneID), "/"),
                  function(g) paste(sort(g), collapse = ","), "")
    padj <- if (!is.null(df$p.adjust)) df$p.adjust else rep(1, nrow(df))
    ord <- order(padj)
    df <- df[ord, , drop = FALSE]
    key <- key[ord]
    ## One row per distinct gene set; the best-corrected name represents it.
    out <- df[!duplicated(key), , drop = FALSE]
    keybest <- key[!duplicated(key)]
    all_names <- split(df$Description, key)
    out$aliases <- vapply(keybest, function(k) {
        first <- out$Description[match(k, keybest)]
        others <- setdiff(all_names[[k]], first)
        if (!length(others)) NA_character_ else
            paste(utils::head(others, 5), collapse = "; ")
    }, "")
    rownames(out) <- NULL
    out
}

## Test every collection against every pathway.
##
## The denominator is the caller's universe, not the pathway's own genes: the
## question is whether this regulator's targets are over-represented *among the
## pathway's genes relative to the whole experiment*, which is the same
## background the collection layer used. Using the pathway as its own background
## would make every set look enriched.
##
## Correction is within the pathway, across the sets tested for it -- the same
## rule padj_within_source follows, applied to the smaller family it belongs to.
.mechanism_pathway_layer <- function(pathways, universe, sets,
                                     minGSSize = 5, maxGSSize = Inf,
                                     top = 20L, top_regulators = 10L) {
    if (is.null(pathways) || !nrow(pathways) || !length(sets)) return(NULL)
    pathways <- utils::head(pathways, max(1L, as.integer(top)))

    rows <- list()
    for (i in seq_len(nrow(pathways))) {
        pgenes <- strsplit(as.character(pathways$geneID[i]), "/")[[1]]
        pgenes <- unique(pgenes[!is.na(pgenes) & nzchar(pgenes)])
        if (length(pgenes) < 2L) next
        for (src in names(sets)) {
            r <- tryCatch(
                suppressWarnings(enricher(
                    pgenes, universe = universe, gson = sets[[src]],
                    minGSSize = minGSSize, maxGSSize = maxGSSize,
                    pvalueCutoff = 1, qvalueCutoff = 1)),
                error = function(e) NULL)
            if (is.null(r)) next
            d <- as.data.frame(r)
            if (!nrow(d)) next
            d <- d[!is.na(d$ID), , drop = FALSE]
            if (!nrow(d)) next
            rows[[length(rows) + 1L]] <- data.frame(
                pathway             = as.character(pathways$ID[i]),
                pathway_description = as.character(pathways$Description[i] %||%
                                                   pathways$ID[i]),
                pathway_direction   = as.character(pathways$direction[i] %||%
                                                   NA_character_),
                n_pathway_genes     = length(pgenes),
                ## The pathway's own enrichment p-value, kept apart from
                ## padj_within_pathway (which is the attribution's). Without it a
                ## reader has no way to tell a strongly enriched pathway that
                ## attracted a weak attribution from the reverse, and the
                ## narrative context cannot order by what the enrichment ranked
                ## highest.
                pathway_padj        = as.numeric(pathways$p.adjust[i] %||% NA_real_),
                regulator           = mechanism_alias_regulator(
                                          mechanism_entity(d$ID, src)),
                source              = src,
                term                = d$ID,
                Count               = as.integer(d$Count),
                setSize             = as.integer(d$SetSize),
                padj_within_pathway = d$p.adjust,
                geneID              = d$geneID,
                stringsAsFactors    = FALSE)
        }
    }
    if (!length(rows)) return(NULL)
    out <- do.call(rbind, rows)

    ## Every sizeable set is "enriched" against a thirty-gene pathway, for the
    ## same arithmetic that makes any big set look significant on a small list:
    ## a 1300-gene set against 14000 background genes has an expected overlap of
    ## three, so almost any overlap clears the test. Reporting all of them gave
    ## fourteen thousand rows for twenty pathways, which is not an answer. Two
    ## cuts make it one: keep the corrected-significant rows, then the strongest
    ## few per pathway.
    ## A consensus sequence names no factor. It is the single most significant
    ## set in most lists, so leaving it in puts "GTGGGTGK_UNKNOWN regulates
    ## blood vessel development" at the top of every pathway -- a sentence that
    ## names nothing and cannot be checked against any perturbation dataset.
    ## The same exclusion the verdict layer already applies.
    unnamed <- mechanism_unnamed(out$regulator)
    out <- out[!unnamed, , drop = FALSE]

    out <- out[!is.na(out$padj_within_pathway) & out$padj_within_pathway < 0.05, ,
               drop = FALSE]
    if (!nrow(out)) return(NULL)
    out <- out[order(out$pathway, out$padj_within_pathway), , drop = FALSE]
    ## One regulator is one claim about one pathway. C3:TFT holds near-duplicate
    ## sets for a factor -- HIF1_Q3 and HIF1_Q5 are the same genes under two
    ## motif names -- so without this the same fact is reported twice and looks
    ## like two independent findings.
    out <- out[!duplicated(paste(out$pathway, out$regulator)), , drop = FALSE]
    per_pathway <- max(1L, as.integer(top_regulators))
    out <- do.call(rbind, lapply(split(out, out$pathway), utils::head, n = per_pathway))
    rownames(out) <- NULL
    out
}

## Pathway attributions as narrative context.
##
## The pathway layer answers "which regulator explains this pathway", which is
## worth writing down -- but it is a re-reading of the same collections the
## verdicts already came from, so it must never add an entity or strengthen a
## claim. Two guards make that structural rather than a request:
##
##   * only regulators the adjudication already admitted are listed, so the
##     context cannot name anything the verdict table did not;
##   * the block says so in its own text, because the instruction has to travel
##     with the data. Putting it in the system prompt instead would change the
##     prompt of the four legacy input types, whose bytes are frozen.
.mechanism_pathway_context <- function(adj, max_pathways = 10L,
                                       max_regulators = 3L) {
    pw <- adj$evidence$pathways
    if (is.null(pw) || !nrow(pw)) return(NULL)

    admitted <- adj$verdicts$regulator[
        adj$verdicts$verdict %in% c("supported", "single_source")]
    if (!length(admitted)) return(NULL)
    pw <- pw[pw$regulator %in% admitted, , drop = FALSE]
    if (!nrow(pw)) return(NULL)

    ## Ordered by the pathway's best corrected p-value, not by its identifier.
    ## The layer returns its rows grouped by pathway id (split() sorts the
    ## levels), so taking the first few in row order would show whichever
    ## identifiers sort early -- for GO, essentially at random -- instead of the
    ## pathways the enrichment actually ranked highest.
    ## Ordered by the pathway's own enrichment p-value where it is known: the
    ## section is about pathways, so the reader expects the ones the enrichment
    ## ranked highest, not the ones that happen to carry the sharpest
    ## attribution. Falls back to the attribution when the pathway p is missing,
    ## which is the case for a caller-supplied result that omits it.
    order_key <- if (!is.null(pw$pathway_padj) && any(!is.na(pw$pathway_padj))) {
        tapply(pw$pathway_padj, pw$pathway, min, na.rm = TRUE)
    } else {
        tapply(pw$padj_within_pathway, pw$pathway, min, na.rm = TRUE)
    }
    ids <- utils::head(names(sort(order_key)), max_pathways)
    lines <- vapply(ids, function(id) {
        rows <- pw[pw$pathway == id, , drop = FALSE]
        rows <- utils::head(rows[order(rows$padj_within_pathway), , drop = FALSE],
                            max_regulators)
        direction <- rows$pathway_direction[1]
        if (is.na(direction)) direction <- "direction not stated"
        sprintf("- %s (%s): %s", rows$pathway_description[1], direction,
                paste0(rows$regulator, " padj ", signif(rows$padj_within_pathway, 2),
                       collapse = "; "))
    }, "")

    c(paste0("These are attributions of entities this analysis has already graded. ",
             "They are context, not evidence: do not introduce a regulator that is ",
             "not already supported or single_source, and do not describe a claim as ",
             "stronger because a pathway is attributed to it."),
      "",
      lines)
}
