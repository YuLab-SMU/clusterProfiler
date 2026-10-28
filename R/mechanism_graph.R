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
#' most results and including them would bury the findings; [adjudicate()] lists
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
#' @seealso [plot_mechanism()], [evidence()], [adjudicate()]
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

    findings <- consensus_findings(adjudicate(x))
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
            i = "`adjudicate()` lists what was examined and why it did not qualify."
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
        from    = members$regulator,
        to      = members$gene,
        type    = "target",
        term    = members$term,
        level   = members$level,
        source  = members$source,
        ## Carried on the edge so attribute_genes() can rank competing claims
        ## without re-joining to the verdict table.
        padj    = members$padj,
        n_genes = members$n_genes,
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
                        source = character(), padj = numeric(),
                        n_genes = integer(), stringsAsFactors = FALSE)
    if (!nrow(findings)) return(empty)

    r <- x$regulators
    r <- r[r$regulator %in% findings$regulator, , drop = FALSE]
    r <- r[mechanism_significant(r), , drop = FALSE]
    sources <- intersect(unique(r$source),
                         c("tf", "tfunion", "mirna", "cellstate", "tissue"))
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
                padj      = r$padj_within_source[match(term, r$term)],
                n_genes   = length(genes),
                stringsAsFactors = FALSE)
        }
    }
    if (!length(out)) return(empty)
    do.call(rbind, out)
}

#' Attribute each gene to the regulator that best explains it
#'
#' @description
#' Collapses the many-to-many regulator/gene relation into a tree: each input
#' gene is assigned to exactly one regulator, and the assignment follows the
#' evidence rather than the layout.
#'
#' @details
#' A network answers "who claims this gene". A mechanism needs "who explains
#' it", and when three regulators each list the same gene the figure cannot say
#' which one the evidence favours -- every line looks the same. This assigns one
#' parent per gene, by a rule stated in order: strongest grade first, then the
#' most significant set, then the largest overlap. The result is a forest that
#' can be drawn as a tree, with the regulators at the root end and the genes at
#' the tips.
#'
#' A gene claimed by no regulator is absent, not unassigned: it is left out
#' rather than attached to a root it has nothing to do with.
#'
#' The order matters and is deliberate. Grade before significance because a
#' measured mechanism outranks a well-powered prediction; significance before
#' overlap because a set that barely reaches threshold should not win on size.
#'
#' @param x A `mechgraph` from [mechanism_graph()].
#' @return A data frame with `gene`, `regulator`, `level`, `padj` and
#'   `n_claimed` (how many regulators listed it, which is the ambiguity this
#'   rule resolves).
#' @seealso [mechanism_graph()], [plot_mechanism_tree()]
#' @examples
#' \dontrun{
#' attribute_genes(mechanism_graph(ev))
#' }
#' @export
attribute_genes <- function(x) {
    if (!inherits(x, "mechgraph")) {
        rlang::abort("`x` must be a mechgraph from mechanism_graph().")
    }
    e <- x$edges
    if (is.null(e) || !nrow(e)) {
        return(data.frame(gene = character(), regulator = character(),
                          level = character(), padj = numeric(),
                          n_claimed = integer(), stringsAsFactors = FALSE))
    }
    nodes <- x$nodes
    e$grade <- nodes$grade[match(e$from, nodes$id)]
    e$padj <- suppressWarnings(as.numeric(e$padj))
    e$n_genes <- suppressWarnings(as.numeric(e$n_genes))

    ## A family label cannot be the answer to "which gene does this", so it
    ## loses to a real regulator at the same grade.
    e$is_family <- nodes$type[match(e$from, nodes$id)] == "family"
    grade_rank <- c(A = 1L, B = 2L, C = 3L)
    e$rank <- grade_rank[e$grade]
    e$rank[is.na(e$rank)] <- 4L

    claimed <- table(e$to)
    e <- e[order(e$to, e$rank, e$is_family,
                 ifelse(is.na(e$padj), Inf, e$padj),
                 -ifelse(is.na(e$n_genes), 0, e$n_genes)), , drop = FALSE]
    e <- e[!duplicated(e$to), , drop = FALSE]

    out <- data.frame(
        gene      = e$to,
        regulator = e$from,
        level     = e$level,
        padj      = e$padj,
        n_claimed = as.integer(claimed[e$to]),
        stringsAsFactors = FALSE)
    out[order(out$regulator, out$gene), , drop = FALSE]
}

#' Draw the mechanism as a tree
#'
#' @description
#' Renders the attribution from [attribute_genes()] with ggtree: one branch per
#' regulator, its assigned genes as tips, coloured by the grade that won the
#' assignment.
#'
#' @details
#' A tree rather than a network because the question a mechanism has to answer
#' is directional -- which factor accounts for which genes -- and a network
#' draws every claim with the same weight, leaving the reader to resolve
#' overlaps by eye. The tree states the resolution: when three regulators list
#' the same gene, one of them is named.
#'
#' The tree is drawn as a fan when there are many regulators and as a rectangle
#' otherwise, because a fan with three branches is mostly empty space.
#'
#' @param x A `mechgraph` from [mechanism_graph()], or the data frame from
#'   [attribute_genes()].
#' @param layout Passed to [ggtree::ggtree()]. The default `"tree_and_leaf"`
#'   puts the regulators in an inner ring and their genes outside, which keeps
#'   the two sets of labels apart; `"rectangular"` reads top to bottom and suits
#'   a few regulators; `"fan"` and `"dendrogram"` are also available.
#' @param ... Passed to [ggtree::ggtree()].
#' @return A `ggtree` plot.
#' @seealso [attribute_genes()], [plot_mechanism()]
#' @examples
#' \dontrun{
#' plot_mechanism_tree(mechanism_graph(ev))
#' }
#' @export
plot_mechanism_tree <- function(x,
                                layout = c("tree_and_leaf", "dendrogram",
                                           "rectangular", "slanted",
                                           "fan", "circular"), ...) {
    layout <- match.arg(layout)
    for (pkg in c("ggtree", "ape")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            rlang::abort(paste0("Package '", pkg, "' is required to draw a tree."))
        }
    }

    assign <- if (inherits(x, "mechgraph")) attribute_genes(x) else x
    if (!is.data.frame(assign) || !nrow(assign)) {
        rlang::abort(c(
            "There is nothing to draw: no gene is attributed to a regulator.",
            i = "`attribute_genes()` returns the assignment this plot is made from."
        ))
    }

    phy <- .attribution_to_phylo(assign)

    ## ggtree's %<+% joins a data frame onto the plot's node table by `label`.
    ## as.treedata(phy, data = ) does not do this -- it accepts the argument and
    ## drops it -- which is why aes(colour = grade) kept failing on a column
    ## that existed in tip_data and nowhere in the plot.
    p <- ggtree::`%<+%`(ggtree::ggtree(phy, layout = layout, ...), phy$tip_data) +
        ggtree::geom_tiplab(size = 2, offset = 0.02) +
        ggtree::geom_tippoint(ggplot2::aes(colour = .data$grade), size = 1.5) +
        ## The regulators sit on internal nodes: without naming them the figure
        ## shows which genes grouped but not who grouped them. node.label is
        ## what ggtree reads for internal nodes, so the names go there rather
        ## than into a separate frame with no coordinates.
        ## The regulator name sits at the node it belongs to. A plain text label
        ## sat on top of the incoming branch; geom_label2 draws a box behind it,
        ## which is what makes the name readable where the branches converge.
        ggtree::geom_label2(
            ggplot2::aes(label = .data$label, subset = !.data$isTip),
            size = 3, colour = "grey10", fill = "white",
            label.size = 0, label.padding = grid::unit(0.06, "lines"),
            nudge_x = -0.02, hjust = 1) +
        ggplot2::scale_colour_manual(
            values = c(A = "#B2182B", B = "#EF8A62", C = "#FDDBC7"),
            name = NULL, drop = FALSE) +
        ## The branches are a drawing device, not a distance: nothing here was
        ## estimated as one, so the axis that would imply it is removed.
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(colour = NULL)
    p
}

## The attribution as a phylo, via a shortest-path tree.
##
## The assigned edges already form a forest -- each gene has exactly one parent
## -- and adding a super-root above the regulators makes it a tree. igraph's
## bfs() then gives the parent of every node, which is the tree. This replaces a
## hand-written Newick string: same result, and it uses the graph library's own
## traversal instead of a string format that has to be re-parsed.
##
## The tip label a reader sees is the gene; which regulator it hangs from is
## said by the branch, and repeating it in every tip would double the width of
## the figure.
.attribution_to_phylo <- function(assign) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        rlang::abort("Package 'igraph' is required to build the tree.")
    }
    assign <- assign[order(assign$regulator, assign$gene), , drop = FALSE]
    ROOT <- "root"

    edges <- rbind(
        data.frame(from = assign$regulator, to = assign$gene,
                   stringsAsFactors = FALSE),
        data.frame(from = ROOT, to = unique(assign$regulator),
                   stringsAsFactors = FALSE))
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    ## bfs() is the shortest-path tree: from the super-root every regulator is
    ## one hop away and every gene two, so no gene can be reached through
    ## another gene.
    tree <- igraph::bfs(g, root = ROOT, parent = TRUE)

    parent <- tree$parent
    ids <- names(parent)
    ## Drop the super-root: ape wants a tree whose root is a real node, and the
    ## root is an artefact of the traversal.
    keep <- ids != ROOT
    ids <- ids[keep]
    parent <- as.character(parent[keep])
    parent[parent == ROOT] <- NA_character_

    ## ape numbers tips first, then internal nodes; build the edge list in that
    ## order so the phylo is valid.
    is_tip <- vapply(ids, function(i) !(i %in% parent), TRUE)
    tips <- ids[is_tip]
    internal <- ids[!is_tip]
    number <- setNames(seq_along(c(tips, internal)), c(tips, internal))

    ## An internal node with a single child is a chain and will not render; the
    ## regulators all have at least one assigned gene, which is what makes this
    ## safe rather than merely likely.
    edge <- do.call(rbind, lapply(ids, function(i) {
        p <- parent[[i]]
        if (is.na(p)) return(NULL)
        data.frame(from = number[[p]], to = number[[i]], stringsAsFactors = FALSE)
    }))
    if (is.null(edge) || !nrow(edge)) {
        rlang::abort("The attribution produced no tree: no gene has a parent.")
    }

    phy <- structure(
        list(edge = as.matrix(edge),
             tip.label = tips,
             Nnode = length(internal),
             edge.length = rep(1, nrow(edge))),
        class = "phylo")
    phy <- .phylo_reorder(phy)

    reg_of <- assign$regulator[match(phy$tip.label, assign$gene)]
    phy$tip_data <- data.frame(
        label = phy$tip.label,
        regulator = reg_of,
        grade = assign$level[match(phy$tip.label, assign$gene)],
        stringsAsFactors = FALSE)

    ## Name the internal node carrying each regulator. The super-root is left
    ## unnamed because it is an artefact of the traversal, not an entity.
    n_tip <- length(phy$tip.label)
    internal <- sort(unique(phy$edge[, 1][phy$edge[, 1] > n_tip]))
    node_label <- rep("", length(internal))
    names(node_label) <- as.character(internal)
    for (reg in unique(assign$regulator)) {
        tips <- which(reg_of == reg)
        if (!length(tips)) next
        parents <- phy$edge[phy$edge[, 2] %in% tips, 1]
        if (!length(parents)) next
        parent <- names(sort(table(parents), decreasing = TRUE))[1]
        node_label[[parent]] <- reg
    }
    phy$node.label <- unname(node_label)
    phy
}

## Name the internal node that carries each regulator.
##
## The tree is unlabelled on the inside by construction, so without this a
## reader sees which genes grouped together but not who grouped them -- which is
## the one thing the figure exists to say. The super-root is left unnamed
## because it is an artefact of the traversal, not an entity.
.mechanism_node_labels <- function(phy, assign) {
    n_tip <- length(phy$tip.label)
    regs <- unique(assign$regulator)
    ## Each regulator is the parent of its own tips, so its node is the one
    ## whose children are all tips of that regulator.
    out <- data.frame(node = integer(), label = character(),
                      stringsAsFactors = FALSE)
    for (reg in regs) {
        tips <- which(assign$regulator[match(phy$tip.label, assign$gene)] == reg)
        if (!length(tips)) next
        parents <- phy$edge[phy$edge[, 2] %in% tips, 1]
        parent <- names(sort(table(parents), decreasing = TRUE))[1]
        out <- rbind(out, data.frame(node = as.integer(parent), label = reg,
                                     stringsAsFactors = FALSE))
    }
    out
}

## Number the internal nodes and order the edges so ape's invariants hold.
.phylo_reorder <- function(phy) {
    n_tip <- length(phy$tip.label)
    internal <- setdiff(phy$edge[, 1], seq_len(n_tip))
    ## Root first, then the rest, which is the order ape expects.
    if (!identical(internal[1], n_tip + 1L)) {
        map <- setNames(c(seq_len(n_tip), n_tip + seq_along(internal)),
                        c(seq_len(n_tip), internal))
        phy$edge[] <- map[as.character(phy$edge)]
    }
    phy
}

#' Attribute each gene to the regulator that best explains it
#'
#' @description
#' Collapses the many-to-many regulator/gene relation into a tree: each input
#' gene is assigned to exactly one regulator, and the assignment follows the
#' evidence rather than the layout.
#'
#' @details
#' A network answers "who claims this gene". A mechanism needs "who explains
#' it", and when three regulators each list the same gene the figure cannot say
#' which one the evidence favours -- every line looks the same. This assigns one
#' parent per gene, by a rule stated in order: strongest grade first, then the
#' most significant set, then the largest overlap. The result is a forest that
#' can be drawn as a tree, with the regulators at the root end and the genes at
#' the tips.
#'
#' A gene claimed by no regulator is absent, not unassigned: it is left out
#' rather than attached to a root it has nothing to do with.
#'
#' The order matters and is deliberate. Grade before significance because a
#' measured mechanism outranks a well-powered prediction; significance before
#' overlap because a set that barely reaches threshold should not win on size.
#'
#' @param x A `mechgraph` from [mechanism_graph()].
#' @return A data frame with `gene`, `regulator`, `level`, `padj` and
#'   `n_claimed` (how many regulators listed it, which is the ambiguity this
#'   rule resolves).
#' @seealso [mechanism_graph()], [plot_mechanism_tree()]
#' @examples
#' \dontrun{
#' attribute_genes(mechanism_graph(ev))
#' }
#' @export
attribute_genes <- function(x) {
    if (!inherits(x, "mechgraph")) {
        rlang::abort("`x` must be a mechgraph from mechanism_graph().")
    }
    e <- x$edges
    if (is.null(e) || !nrow(e)) {
        return(data.frame(gene = character(), regulator = character(),
                          level = character(), padj = numeric(),
                          n_claimed = integer(), stringsAsFactors = FALSE))
    }
    nodes <- x$nodes
    e$grade <- nodes$grade[match(e$from, nodes$id)]
    e$padj <- suppressWarnings(as.numeric(e$padj))
    e$n_genes <- suppressWarnings(as.numeric(e$n_genes))

    ## A family label cannot be the answer to "which gene does this", so it
    ## loses to a real regulator at the same grade.
    e$is_family <- nodes$type[match(e$from, nodes$id)] == "family"
    grade_rank <- c(A = 1L, B = 2L, C = 3L)
    e$rank <- grade_rank[e$grade]
    e$rank[is.na(e$rank)] <- 4L

    claimed <- table(e$to)
    e <- e[order(e$to, e$rank, e$is_family,
                 ifelse(is.na(e$padj), Inf, e$padj),
                 -ifelse(is.na(e$n_genes), 0, e$n_genes)), , drop = FALSE]
    e <- e[!duplicated(e$to), , drop = FALSE]

    out <- data.frame(
        gene      = e$to,
        regulator = e$from,
        level     = e$level,
        padj      = e$padj,
        n_claimed = as.integer(claimed[e$to]),
        stringsAsFactors = FALSE)
    out[order(out$regulator, out$gene), , drop = FALSE]
}

#' Draw the mechanism as a tree
#'
#' @description
#' Renders the attribution from [attribute_genes()] with ggtree: one branch per
#' regulator, its assigned genes as tips, coloured by the grade that won the
#' assignment.
#'
#' @details
#' A tree rather than a network because the question a mechanism has to answer
#' is directional -- which factor accounts for which genes -- and a network
#' draws every claim with the same weight, leaving the reader to resolve
#' overlaps by eye. The tree states the resolution: when three regulators list
#' the same gene, one of them is named.
#'
#' The tree is drawn as a fan when there are many regulators and as a rectangle
#' otherwise, because a fan with three branches is mostly empty space.
#'
#' @param x A `mechgraph` from [mechanism_graph()], or the data frame from
#'   [attribute_genes()].
#' @param layout Passed to [ggtree::ggtree()]. The default `"tree_and_leaf"`
#'   puts the regulators in an inner ring and their genes outside, which keeps
#'   the two sets of labels apart; `"rectangular"` reads top to bottom and suits
#'   a few regulators; `"fan"` and `"dendrogram"` are also available.
#' @param ... Passed to [ggtree::ggtree()].
#' @return A `ggtree` plot.
#' @seealso [attribute_genes()], [plot_mechanism()]
#' @examples
#' \dontrun{
#' plot_mechanism_tree(mechanism_graph(ev))
#' }
#' @export
plot_mechanism_tree <- function(x,
                                layout = c("tree_and_leaf", "dendrogram",
                                           "rectangular", "slanted",
                                           "fan", "circular"), ...) {
    layout <- match.arg(layout)
    for (pkg in c("ggtree", "ape")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            rlang::abort(paste0("Package '", pkg, "' is required to draw a tree."))
        }
    }

    assign <- if (inherits(x, "mechgraph")) attribute_genes(x) else x
    if (!is.data.frame(assign) || !nrow(assign)) {
        rlang::abort(c(
            "There is nothing to draw: no gene is attributed to a regulator.",
            i = "`attribute_genes()` returns the assignment this plot is made from."
        ))
    }

    phy <- .attribution_to_phylo(assign)

    ## ggtree's %<+% joins a data frame onto the plot's node table by `label`.
    ## as.treedata(phy, data = ) does not do this -- it accepts the argument and
    ## drops it -- which is why aes(colour = grade) kept failing on a column
    ## that existed in tip_data and nowhere in the plot.
    p <- ggtree::`%<+%`(ggtree::ggtree(phy, layout = layout, ...), phy$tip_data) +
        ggtree::geom_tiplab(size = 2, offset = 0.02) +
        ggtree::geom_tippoint(ggplot2::aes(colour = .data$grade), size = 1.5) +
        ## The regulators sit on internal nodes: without naming them the figure
        ## shows which genes grouped but not who grouped them. node.label is
        ## what ggtree reads for internal nodes, so the names go there rather
        ## than into a separate frame with no coordinates.
        ## The regulator name sits at the node it belongs to. A plain text label
        ## sat on top of the incoming branch; geom_label2 draws a box behind it,
        ## which is what makes the name readable where the branches converge.
        ggtree::geom_label2(
            ggplot2::aes(label = .data$label, subset = !.data$isTip),
            size = 3, colour = "grey10", fill = "white",
            label.size = 0, label.padding = grid::unit(0.06, "lines"),
            nudge_x = -0.02, hjust = 1) +
        ggplot2::scale_colour_manual(
            values = c(A = "#B2182B", B = "#EF8A62", C = "#FDDBC7"),
            name = NULL, drop = FALSE) +
        ## The branches are a drawing device, not a distance: nothing here was
        ## estimated as one, so the axis that would imply it is removed.
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(colour = NULL)
    p
}

## A phylo whose tips are "regulator|gene", so each gene sits under the
## regulator it was attributed to.
##
## Built as a Newick string and parsed by ape rather than assembled from a
## parent vector: this code path is exercised once per figure, and the string
## form is what can be read back when a figure looks wrong.
.attribution_to_phylo <- function(assign) {
    assign <- assign[order(assign$regulator, assign$gene), , drop = FALSE]
    by_reg <- split(assign, assign$regulator, drop = TRUE)

    branches <- vapply(names(by_reg), function(reg) {
        genes <- by_reg[[reg]]$gene
        tips <- paste0("'", reg, "|", genes, "':1")
        paste0("(", paste(tips, collapse = ","), ")")
    }, character(1))

    nw <- paste0("(", paste(branches, collapse = ","), ");")
    phy <- ape::read.tree(text = nw)
    ## Edge lengths are a drawing device, not a measurement: nothing here was
    ## estimated as a distance, so they must not be read as one.
    phy$edge.length <- rep(1, length(phy$edge.length))

    tip_reg <- gsub("'", "", sub("^'?([^|]+)\\|.*$", "\\1", phy$tip.label))
    tip_gene <- gsub("'", "", sub("^.*\\|'?(.*?)'?$", "\\1", phy$tip.label))

    ## The tip label a reader sees is the gene; which regulator it sits under is
    ## already said by the branch it hangs from, so repeating the regulator in
    ## every tip would double the width of the figure.
    phy$tip.label <- tip_gene
    phy$tip_data <- data.frame(
        label = phy$tip.label,
        regulator = tip_reg,
        grade = assign$level[match(tip_gene, assign$gene)],
        stringsAsFactors = FALSE)

    ## Name the internal node carrying each regulator. The super-root is left
    ## unnamed because it is an artefact of the traversal, not an entity.
    n_tip <- length(phy$tip.label)
    internal <- sort(unique(phy$edge[, 1][phy$edge[, 1] > n_tip]))
    node_label <- rep("", length(internal))
    names(node_label) <- as.character(internal)
    for (reg in unique(assign$regulator)) {
        tips <- which(tip_reg == reg)
        if (!length(tips)) next
        parents <- phy$edge[phy$edge[, 2] %in% tips, 1]
        if (!length(parents)) next
        parent <- names(sort(table(parents), decreasing = TRUE))[1]
        if (parent %in% names(node_label)) node_label[[parent]] <- reg
    }
    phy$node.label <- unname(node_label)
    phy
}
