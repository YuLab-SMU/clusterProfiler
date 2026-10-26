#' The mechanism the evidence permits, as a mechgraph
#'
#' @description
#' Turns an `evidenceBundle` into a `mechgraph` object: candidate regulators on
#' one side, the input genes they are credited with on the other, every node and
#' edge carrying the grade it earned. Draw it with `plot()` (rendered by
#' ggtangle) or keep it as data.
#'
#' @details
#' A `mechgraph` rather than a bare data frame, because the object already
#' carries the contract this analysis needs: a node has an `id`, a `type` and a
#' `label`; an edge has a `from`, a `to` and a `type`. Anything else rides along
#' as a column and is passed through to igraph, so the grade travels with the
#' node instead of being re-joined at plot time.
#'
#' A link means **"this collection lists this gene"**, not "this factor
#' regulates this gene". The collection is a target repertoire, and which of its
#' members are real targets is exactly what the analysis could not settle. That
#' distinction is why an edge here is evidence of association regardless of the
#' regulator's own grade.
#'
#' Absent and untested entities are not included. They are the larger part of
#' most results and including them would bury the findings; [consensus()] lists
#' them.
#'
#' @param x An `evidenceBundle` from [evidence()].
#' @param max_genes Most genes to keep per regulator, by overlap size. A
#'   collection term can list a thousand genes and the figure stops being
#'   readable long before that.
#' @param species Species whose collections supply the membership.
#' @param release MSigDB release, taken from the bundle when `NULL`.
#' @return A `mechgraph` object. `nodes` carry `id`, `type` (`"regulator"`,
#'   `"family"` or `"gene"`), `label`, `grade` and `verdict`; `edges` carry
#'   `from`, `to`, `type` (`"target"`), `term`, `level` and `source`.
#'   `metadata` records the input size and the species.
#' @seealso [plot_mechanism()], [evidence()], [consensus()]
#' @examples
#' \dontrun{
#' ev <- evidence(genes, upstream = c("tf", "perturbagen"),
#'                design = list(treatment = "knockdown"))
#' mg <- mechanism_graph(ev)
#' plot(mg)
#' }
#' @export
mechanism_graph <- function(x, max_genes = 8, species = NULL, release = NULL) {
    if (!requireNamespace("mechgraph", quietly = TRUE)) {
        rlang::abort(c(
            "Package 'mechgraph' is required to build a mechanism graph.",
            i = "Install it with remotes::install_github(\"YuLab-SMU/mechgraph\")."
        ))
    }
    if (!inherits(x, "evidenceBundle")) {
        rlang::abort("`x` must be an evidenceBundle from evidence().")
    }

    findings <- consensus_findings(consensus(x))
    if (is.null(species)) species <- x$input$design$species %||% "Homo sapiens"
    if (is.null(release)) {
        rel <- unique(stats::na.omit(x$provenance$release))
        release <- if (length(rel)) rel[1] else "2026.2.Hs"
    }

    ## Two different failures, and reporting one as the other is the mistake
    ## this package keeps having to fix: no finding at all is not the same as
    ## findings whose collections could not be loaded.
    if (!nrow(findings)) {
        rlang::abort(c(
            "There is nothing to draw: no entity reached a supported or single-source verdict.",
            i = "`consensus()` lists what was examined and why it did not qualify."
        ))
    }

    members <- .mechanism_members(x, findings, species, release)
    if (!nrow(members)) {
        rlang::abort(c(
            "The findings are there, but none of them lists a gene from your input.",
            i = "A link is drawn only where a collection contains an input gene.",
            i = "Check `x$diagnostics$failures` for a collection that could not be loaded."))
    }

    ## Cap per regulator, strongest overlap first: a term listing a thousand
    ## genes makes the figure unreadable and hides the finding.
    members <- members[order(members$regulator, -members$n_genes), , drop = FALSE]
    members <- do.call(rbind, lapply(split(members, members$regulator),
                                     utils::head, n = max_genes))
    rownames(members) <- NULL

    regs <- findings[findings$regulator %in% members$regulator, , drop = FALSE]
    regs <- regs[!duplicated(regs$regulator), , drop = FALSE]

    regulator_nodes <- data.frame(
        id      = regs$regulator,
        type    = ifelse(regs$family, "family", "regulator"),
        label   = regs$regulator,
        grade   = regs$max_level,
        verdict = regs$verdict,
        ## A family label is not one gene, so it is a separate type rather than
        ## a regulator of unknown kind.
        evidence_kind = regs$evidence_kind %||% NA_character_,
        stringsAsFactors = FALSE)

    genes <- unique(members$gene)
    gene_nodes <- data.frame(
        id      = genes,
        type    = "gene",
        label   = genes,
        grade   = NA_character_,
        verdict = NA_character_,
        evidence_kind = NA_character_,
        stringsAsFactors = FALSE)

    nodes <- rbind(regulator_nodes, gene_nodes)
    nodes <- nodes[!duplicated(nodes$id), , drop = FALSE]

    edges <- data.frame(
        from   = members$regulator,
        to     = members$gene,
        type   = "target",
        term   = members$term,
        level  = members$level,
        source = members$source,
        stringsAsFactors = FALSE)
    edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, , drop = FALSE]

    mechgraph::mg_from_edges(
        edges = edges,
        nodes = nodes,
        metadata = list(
            source   = "clusterProfiler::evidence",
            n_input  = length(x$input$genes),
            species  = species,
            release  = release,
            ## The caption a reader needs, kept with the data so it cannot be
            ## lost when the figure is redrawn elsewhere.
            note     = paste("Links mean a collection lists the gene, not that",
                             "the factor regulates it.")))
}

#' Draw a mechanism graph
#'
#' @description
#' Renders the `mechgraph` from [mechanism_graph()] with ggtangle, colouring
#' nodes by the grade they earned.
#'
#' @param x A `mechgraph` from [mechanism_graph()].
#' @param layout Layout passed to [ggtangle::ggplot.mechgraph]. `"nicely"` suits
#'   a few findings; `"kk"` spreads a crowded one.
#' @param ... Passed to [ggtangle::ggplot.mechgraph].
#' @return A `ggplot` object.
#' @seealso [mechanism_graph()]
#' @examples
#' \dontrun{
#' plot_mechanism(mechanism_graph(ev))
#' }
#' @export
plot_mechanism <- function(x, layout = "nicely", ...) {
    if (!inherits(x, "mechgraph")) {
        rlang::abort(c(
            "`x` must be a mechgraph from mechanism_graph().",
            i = "Build one with mechanism_graph(evidence(...))."
        ))
    }
    for (pkg in c("ggtangle", "ggplot2")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            rlang::abort(paste0("Package '", pkg, "' is required to draw."))
        }
    }

    ## ggtangle::ggplot() prepares the layout and nothing else -- its layers are
    ## added by geom_edge() and geom_point(), which is what the package expects
    ## a caller to do. Calling only the first produced a blank panel with a
    ## caption, which is how this was found.
    ## Nodes are coloured by the grade they earned -- the whole point of the
    ## figure is that a reader can see which findings are measured and which are
    ## predicted without consulting a table. Labels name the regulators only: a
    ## gene label on 68 nodes makes the panel unreadable, and the regulators are
    ## what the figure is about.
    ## ggtangle reads vertex attributes, so the drawing vocabulary is attached
    ## to the graph rather than kept in a parallel frame.
    ig <- igraph::as.igraph(x)
    igraph::vertex_attr(ig, "drawn_as") <- .mechanism_drawn_as(x$nodes)
    igraph::vertex_attr(ig, "short_label") <- .mechanism_labels(x$nodes)
    x <- ig

    ggtangle::ggplot(x, layout = layout, ...) +
        ggtangle::geom_edge(ggplot2::aes(colour = .data$level), alpha = 0.35) +
        ggtangle::geom_point(ggplot2::aes(fill = .data$drawn_as, size = .data$type),
                             shape = 21, colour = "grey25") +
        ggtangle::geom_label(ggplot2::aes(label = .data$short_label),
                             size = 2.6, fill = "white", alpha = 0.75,
                             label.size = 0, label.padding = grid::unit(0.08, "lines"),
                             check_overlap = TRUE) +
        ggplot2::scale_size_manual(values = c(regulator = 4.5, family = 4.5, gene = 1.8),
                                   guide = "none") +
        ggplot2::scale_fill_manual(
            values = c("level A" = "#B2182B", "level B" = "#EF8A62",
                       "level C" = "#FDDBC7", "family / complex" = "#7F7F7F",
                       "gene" = "#9ECAE1"),
            name = NULL, drop = FALSE) +
        ggplot2::labs(caption = .mechanism_caption(x), colour = NULL) +
        ggplot2::theme(legend.position = "bottom")
}

## What a node is drawn as, in the vocabulary a reader needs.
##
## The grade decides the colour; a family label is its own class because it is
## not one gene and cannot be one, and a gene carries no grade of its own.
.mechanism_drawn_as <- function(nodes) {
    out <- rep("gene", nrow(nodes))
    out[nodes$type == "family"] <- "family / complex"
    reg <- nodes$type == "regulator"
    grade <- nodes$grade
    out[reg & !is.na(grade)] <- paste0("level ", grade[reg & !is.na(grade)])
    out[reg & is.na(grade)] <- "level C"
    out
}

## Label a regulator, and leave a gene unlabelled.
.mechanism_labels <- function(nodes) {
    ifelse(nodes$type == "gene", "", nodes$label)
}

.mechanism_caption <- function(x) {
    note <- x$metadata$note %||%
        "Links mean a collection lists the gene, not that the factor regulates it."
    paste(strwrap(note, width = 90), collapse = "\n")
}

## Regulator -> gene links, taken from the collection that produced the finding.
##
## Only genes the caller supplied are kept: the question is what of *this* list
## each regulator explains, not what else the collection contains.
.mechanism_members <- function(x, findings, species, release) {
    empty <- data.frame(regulator = character(), gene = character(),
                        term = character(), level = character(),
                        source = character(), n_genes = integer(),
                        stringsAsFactors = FALSE)
    if (!nrow(findings)) return(empty)

    r <- x$regulators
    r <- r[r$regulator %in% findings$regulator, , drop = FALSE]
    r <- r[mechanism_significant(r), , drop = FALSE]
    sources <- intersect(unique(r$source), c("tf", "mirna", "cellstate", "tissue"))
    if (!length(sources) || !nrow(r)) return(empty)

    selected <- unique(as.character(x$input$genes))
    out <- list()
    for (src in sources) {
        g <- tryCatch(mechanism_sets(src, species = species, release = release),
                      error = function(e) NULL)
        if (is.null(g)) next
        m <- g@gsid2gene
        for (term in unique(r$term[r$source == src])) {
            genes <- intersect(unique(m$gene[m$gsid == term]), selected)
            if (!length(genes)) next
            out[[length(out) + 1L]] <- data.frame(
                regulator = rep(mechanism_entity(term, src, species), length(genes)),
                gene      = genes,
                term      = term,
                level     = mechanism_set_level(src, term),
                source    = src,
                n_genes   = length(genes),
                stringsAsFactors = FALSE)
        }
    }
    if (!length(out)) return(empty)
    do.call(rbind, out)
}
