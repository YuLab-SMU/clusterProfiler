##' plot GO graph
##'
##'
##' @title plotGOgraph
##' @param x output of enrichGO or gseGO
##' @param firstSigNodes number of significant nodes (retangle nodes in the graph)
##' @param useInfo additional info
##' @param sigForAll if TRUE the score/p-value of all nodes in the DAG is shown, otherwise only score will be shown
##' @param useFullNames logical
##' @param ... additional parameter of showSigOfNodes, please refer to topGO
##' @return GO DAG graph
## @importClassesFrom topGO topGOdata
## @importFrom topGO showSigOfNodes
## @importFrom topGO annFUN.gene2GO
## @importFrom topGO groupGOTerms
##' @export
##' @seealso
##' \link[topGO]{showSigOfNodes}
##' @author Guangchuang Yu
plotGOgraph <- function(x,
                        firstSigNodes=10,
                        useInfo="all",
                        sigForAll=TRUE,
                        useFullNames=TRUE, ...) {

    ## requireNamespace("topGO") || stop("package topGO is required")
    groupGOTerms <- get_fun_from_pkg("topGO", "groupGOTerms")
    annFUN.gene2GO <- get_fun_from_pkg("topGO", "annFUN.gene2GO")
    showSigOfNodes <- get_fun_from_pkg("topGO", "showSigOfNodes")

    if (! class(x) %in% c("gseaResult", "enrichResult")) {
        stop("x should be output of gseGO or enrichGO...")
    }

    gs <- x@geneSets
    gs.df <- data.frame(gene = unlist(gs),
                        go   = rep(names(gs),
                                   times=sapply(gs, length)))
    gene2GO <- split(as.character(gs.df$go), gs.df$gene)

    if (is(x, "gseaResult")) {
        ont <- x@setType
        allgenes <- x@geneList
        core_genes <- unique(unlist(geneInCategory(x)))
        allgenes[!names(allgenes) %in% core_genes] <- -1
        allgenes[core_genes] <- 1
    } else {
        ont <- x@ontology
        universe <- x@universe
        allgenes <- numeric(length(universe))
        names(allgenes) <- universe
        allgenes[x@gene] <- 1
    }

    selector <- function(scores) return(scores == 1)

    if ( ! ont %in% c("BP", "MF", "CC")) {
        stop("ontology should be one of 'BP', 'MF' or 'CC'...")
    }


    pvalue <- x@result$p.adjust
    names(pvalue) <- x@result$ID

    groupGOTerms()

    GOdata <- new("topGOdata",
                  description="clusterProfiler enrichment results",
                  ontology = ont,
                  allGenes = allgenes,
                  geneSel = selector,
                  annot = annFUN.gene2GO,
                  gene2GO=gene2GO)

    firstSigNodes <- min(firstSigNodes, nrow(x))
    showSigOfNodes(GOdata        = GOdata,
                   termsP.value  = pvalue,
                   firstSigNodes = firstSigNodes,
                   useInfo       = useInfo,
                   sigForAll     = sigForAll,
                   useFullNames  = useFullNames,
                   ...)
}

