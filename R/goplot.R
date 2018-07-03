##' @rdname goplot
##' @exportMethod goplot
setMethod("goplot", signature(x = "enrichResult"),
          function(x, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom="text", ...) {
              goplot.enrichResult(x, showCategory = showCategory,
                                  color = color, layout = layout, geom = geom, ...)
          })

##' @rdname goplot
##' @exportMethod goplot
setMethod("goplot", signature(x = "gseaResult"),
          function(x, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom="text", ...) {
              goplot.enrichResult(x, showCategory = showCategory,
                                  color = color, layout = layout, geom = geom, ...)
          })


##' @rdname goplot
##' @importFrom utils data
##' @import GOSemSim
##' @importFrom ggplot2 scale_fill_gradientn
##' @importFrom grid arrow
##' @importFrom grid unit
##' @import ggraph
##' @importFrom ggraph circle
##' @importFrom ggraph geom_node_label
##' @importFrom AnnotationDbi mget
##' @author guangchuang yu
goplot.enrichResult <- function(x, showCategory = 10, color = "p.adjust", layout = "sugiyama", geom = "text", ...) {
    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    y <- y[1:n,]

    id <- y$ID[1:n]

    if (!exists(".GOSemSimEnv")) GOSemSim_initial()
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    gotbl <- get("gotbl", envir=.GOSemSimEnv)

    GOANCESTOR <- getAncestors(x@ontology)
    anc <- AnnotationDbi::mget(id, GOANCESTOR)
    ca <- anc[[1]]
    for (i in 2:length(anc)) {
        ca <- intersect(ca, anc[[i]])
    }

    uanc <- unique(unlist(anc))
    uanc <- uanc[!uanc %in% ca]
    dag <- gotbl[gotbl$go_id %in% unique(c(id, uanc)),]


    edge <- dag[, c(5, 1, 4)]
    node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
    node$color <- x[node$go_id, color]
    node$size <- sapply(geneSets[node$go_id], length)

    g <- graph.data.frame(edge, directed=TRUE, vertices=node)
    E(g)$relationship <- edge[,3]

    p <- ggraph(g, layout=layout) +
        ## geom_edge_link(aes_(color = ~relationship), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm')) +
        geom_edge_link(aes_(linetype = ~relationship), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm'), colour="darkgrey") +
        ## geom_node_point(size = 5, aes_(fill=~color), shape=21) +
        geom_node_point(size = 5, aes_(color=~color)) +
        theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE))
    ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE))

    if (geom == "label") {
        p <- p + geom_node_label(aes_(label=~Term, fill=~color), repel=TRUE) +
            scale_fill_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE), na.value="white")
        ## scale_fill_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE), na.value='white')
    } else {
        p <- p + geom_node_text(aes_(label=~Term), repel=TRUE)
    }
    return(p)
}


GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
getAncestors <- getFromNamespace("getAncestors", "GOSemSim")
