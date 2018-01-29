##' Compare gene clusters functional profile
##'
##' Given a list of gene set, this function will compute profiles of each gene
##' cluster.
##'
##'
##' @param geneClusters a list of entrez gene id. Alternatively, a formula of type Entrez~group
##' @param fun One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .
##' @param data if geneClusters is a formula, the data from which the clusters must be extracted.
##' @param ...  Other arguments.
##' @return A \code{clusterProfResult} instance.
##' @importFrom methods new
##' @importFrom stats formula
##' @importFrom plyr llply
##' @importFrom plyr ldply
##' @importFrom plyr dlply
##' @importFrom plyr rename
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
##'   \code{\link{enrichGO}}
##' @keywords manip
##' @examples
##' \dontrun{
##' data(gcSample)
##' xx <- compareCluster(gcSample, fun="enrichKEGG",
##'                      organism="hsa", pvalueCutoff=0.05)
##' as.data.frame(xx)
##' # plot(xx, type="dot", caption="KEGG Enrichment Comparison")
##'
##' ## formula interface
##' mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
##'                             '100127206', '100128071'),
##'                    group = c('A', 'A', 'A', 'B', 'B', 'B'),
##'                    othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
##' xx.formula <- compareCluster(Entrez~group, data=mydf,
##'                              fun='groupGO', OrgDb='org.Hs.eg.db')
##' as.data.frame(xx.formula)
##'
##' ## formula interface with more than one grouping variable
##' xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf,
##'                                        fun='groupGO', OrgDb='org.Hs.eg.db')
##' as.data.frame(xx.formula.twogroups)
##' }
compareCluster <- function(geneClusters, fun="enrichGO", data='', ...) {
    fun_name <- fun
    fun <- eval(parse(text=fun))

    # Use formula interface for compareCluster
    if (typeof(geneClusters) == 'language') {
        if (!is.data.frame(data)) {
            stop ('no data provided with formula for compareCluster')
        } else {
            genes.var       = all.vars(geneClusters)[1]
            grouping.formula = gsub('^.*~', '~', as.character(as.expression(geneClusters)))   # For formulas like x~y+z
            geneClusters = dlply(.data=data, formula(grouping.formula), .fun=function(x) {as.character(x[[genes.var]])})
        }
    }
    clProf <- llply(geneClusters,
                    .fun=function(i) {
                        x=suppressMessages(fun(i, ...))
                        if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
                            as.data.frame(x)
                        }
                    }
                    )
    clusters.levels = names(geneClusters)
    clProf.df <- ldply(clProf, rbind)

    if (nrow(clProf.df) == 0) {
        stop("No enrichment found in any of gene cluster, please check your input...")
    }

    clProf.df <- rename(clProf.df, c(.id="Cluster"))
    clProf.df$Cluster = factor(clProf.df$Cluster, levels=clusters.levels)

    if (is.data.frame(data) && grepl('+', grouping.formula)) {
        groupVarName <- strsplit(grouping.formula, split="\\+") %>% unlist %>%
            gsub("~", "", .) %>% gsub("^\\s*", "", .) %>% gsub("\\s*$", "", .)
        groupVars <- sapply(as.character(clProf.df$Cluster), strsplit, split="\\.") %>% do.call(rbind, .)
        for (i in seq_along(groupVarName)) {
            clProf.df[, groupVarName[i]] <- groupVars[,i]
        }
        i <- which(colnames(clProf.df) %in% groupVarName)
        j <- (1:ncol(clProf.df))[-c(1, i)]
        clProf.df <- clProf.df[, c(1, i, j)]
    }

    ##colnames(clProf.df)[1] <- "Cluster"
    new("compareClusterResult",
        compareClusterResult = clProf.df,
        geneClusters = geneClusters,
        fun = fun_name,
        .call = match.call(expand.dots=TRUE)
	)
}



## show method for \code{compareClusterResult} instance
##
##
## @name show
## @alias show
## @docType methods
## @rdname show-methods
##
## @title show method
## @param object A \code{compareClusterResult} instance.
## @return message
## @importFrom methods show
## @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @importFrom utils str
setMethod("show", signature(object="compareClusterResult"),
          function (object){
              cmsg <- paste("  Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.",
                            "  clusterProfiler: an R package for comparing biological themes among",
                            "  gene clusters. OMICS: A Journal of Integrative Biology 2012,",
                            "  16(5):284-287",
                            sep="\n", collapse="\n")

              geneClusterLen <- length(object@geneClusters)
              fun <- object@fun
              result <- object@compareClusterResult
              clusts <- split(result, result$Cluster)
              nterms <- sapply(clusts, nrow)

              cat("#\n# Result of Comparing", geneClusterLen, "gene clusters", "\n#\n")
              cat("#.. @fun", "\t", fun, "\n")
              cat("#.. @geneClusters", "\t")
              str(object@geneClusters)
              cat("#...Result", "\t")
              str(result)
              cat("#.. number of enriched terms found for each gene cluster:\n")
              for (i in seq_along(clusts)) {
                  cat("#..  ", paste0(names(nterms)[i], ":"), nterms[i], "\n")
              }
              cat("#\n#...Citation\n")
              citation_msg <- NULL
              if (fun == "enrichDO" || fun == "enrichNCG") {
                  citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
                                        "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
                                        "  analysis. Bioinformatics 2015 31(4):608-609",
                                        sep="\n", collapse="\n")
              } else if (fun == "enrichPathway") {
                  citation_msg <- paste("  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
                                        "  reactome pathway analysis and visualization. Molecular BioSystems",
                                        "  2016, 12(2):477-479", sep="\n", collapse="\n")
              }
              if (!is.null(citation_msg)) {
                  cat(paste0("1.", citation_msg), "\n\n")
                  cat(paste0("2.", cmsg), "\n\n")
              } else {
                  cat(cmsg, "\n\n")
              }
          })

## summary method for \code{compareClusterResult} instance
##
##
## @name summary
## @alias summary
## @docType methods
## @rdname summary-methods
##
## @title summary method
## @param object A \code{compareClusterResult} instance.
## @return A data frame
## @importFrom stats4 summary
## @exportMethod summary
## @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("summary", signature(object="compareClusterResult"),
          function(object, ...) {
              warning("summary method to convert the object to data.frame is deprecated, please use as.data.frame instead.")
              return(as.data.frame(object, ...))
          }
          )



## ##' @rdname plot-methods
## ##' @importFrom stats4 plot
## ##' @importFrom enrichplot dotplot
## ##' @aliases plot,compareClusterResult,ANY-method
## ##' @param x compareClusterResult object
## ##' @param type one of bar or dot
## ##' @param colorBy one of pvalue or p.adjust
## ##' @param showCategory category numbers
## ##' @param by one of geneRatio, Percentage or count
## ##' @param split ONTOLOGY or NULL
## ##' @param includeAll logical
## ##' @param font.size font size
## ##' @param title figure title
## setMethod("plot", signature(x="compareClusterResult"),
##           function(x,
##                    type="dot",
##                    colorBy="p.adjust",
##                    showCategory=5,
##                    by="geneRatio",
##                    split=NULL,
##                    includeAll=TRUE,
##                    font.size=12,
##                    title=""
##                    ) {
##               if (type == "dot" || type == "dotplot") {
##                   dotplot(x,
##                           colorBy      = colorBy,
##                           showCategory = showCategory,
##                           by           = by,
##                           split        = split,
##                           includeAll   = includeAll,
##                           font.size    = font.size,
##                           title        = title
##                           )
##               } else if (type == "bar" || type == "barplot") {
##                   barplot.compareClusterResult(x, colorBy, showCategory, by, split = split, includeAll, font.size, title)
##               } else {
##                   stop("type should be one of 'dot' or 'bar'...")
##               }
##           })


##' dot plot method
##'
##'
##' @docType methods
##' @title dotplot
##' @rdname dotplot-methods
##' @aliases dotplot,compareClusterResult,ANY-method
##' @param object compareClusterResult object
##' @param x x variable
##' @param color one of pvalue or p.adjust
##' @param showCategory category numbers
##' @param by one of geneRatio, Percentage or count
##' @param split ONTOLOGY or NULL
##' @param includeAll logical
##' @param font.size font size
##' @param title figure title
##' @importFrom enrichplot dotplot
##' @exportMethod dotplot
setMethod("dotplot", signature(object="compareClusterResult"),
          function(object,
                   x = ~Cluster,
                   color ="p.adjust",
                   showCategory=5,
                   split=NULL,
                   font.size=12,
                   title="",
                   by="geneRatio",
                   includeAll=TRUE
                   ) {
              dotplot.compareClusterResult(object, x=x, colorBy = color,
                                           showCategory = showCategory, by = by,
                                           includeAll = includeAll,
                                           split=split, font.size = font.size,
                                           title = title)
          })


barplot.compareClusterResult <- function(height, color="p.adjust", showCategory=5,
                                         by="geneRatio", includeAll=TRUE, font.size=12, title="", ...) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an compareClusterResult object.
    df <- fortify(height, showCategory=showCategory, by=by, includeAll=includeAll)
    plotting.clusterProfile(df, type="bar", colorBy=color, by=by, title=title, font.size=font.size)
}


##' merge a list of enrichResult objects to compareClusterResult
##'
##'
##' @title merge_result
##' @param enrichResultList a list of enrichResult objects
##' @return a compareClusterResult instance
##' @author Guangchuang Yu
##' @importFrom plyr ldply
##' @export
merge_result <- function(enrichResultList) {
    if ( !is(enrichResultList, "list")) {
        stop("input should be a name list...")
    }
    if ( is.null(names(enrichResultList))) {
        stop("input should be a name list...")
    }
    x <- lapply(enrichResultList, as.data.frame)
    names(x) <- names(enrichResultList)
    y <- ldply(x, "rbind")
    y <- rename(y, c(.id="Cluster"))
    y$Cluster = factor(y$Cluster, levels=names(enrichResultList))
    new("compareClusterResult",
        compareClusterResult = y)

}
