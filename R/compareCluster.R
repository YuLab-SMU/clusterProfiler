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
##' @importFrom plyr llply
##' @importFrom plyr ldply
##' @importFrom plyr rename
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
##'   \code{\link{enrichGO}}
##' @keywords manip
##' @examples
##'
##' data(gcSample)
##' xx <- compareCluster(gcSample, fun="enrichKEGG",
##'                      organism="human", pvalueCutoff=0.05)
##' summary(xx)
##' # plot(xx, type="dot", caption="KEGG Enrichment Comparison")
##'
##' ## formula interface
##' mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
##'                             '100127206', '100128071'),
##'                    group = c('A', 'A', 'A', 'B', 'B', 'B'),
##'                    othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
##' xx.formula <- compareCluster(Entrez~group, data=mydf, fun='groupGO')
##' summary(xx.formula)
##'
##' ## formula interface with more than one grouping variable
##' xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf, fun='groupGO')
##' summary(xx.formula.twogroups)
compareCluster <- function(geneClusters, fun="enrichGO", data='', ...) {

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
            x=fun(i, ...)
                        if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
                            summary(x)
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

    ##colnames(clProf.df)[1] <- "Cluster"
    new("compareClusterResult",
        compareClusterResult = clProf.df,
        geneClusters = geneClusters,
        fun = fun
	)
}


##' Class "compareClusterResult"
##' This class represents the comparison result of gene clusters by GO
##' categories at specific level or GO enrichment analysis.
##'
##'
##' @name compareClusterResult-class
##' @aliases compareClusterResult-class show,compareClusterResult-method
##'   summary,compareClusterResult-method plot,compareClusterResult-method
##' @docType class
##' @slot compareClusterResult cluster comparing result
##' @slot geneClusters a list of genes
##' @slot fun one of groupGO, enrichGO and enrichKEGG
##' @exportClass compareClusterResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @exportClass compareClusterResult
##' @seealso \code{\linkS4class{groupGOResult}}
##'   \code{\linkS4class{enrichResult}} \code{\link{compareCluster}}
##' @keywords classes
setClass("compareClusterResult",
         representation = representation(
         compareClusterResult = "data.frame",
         geneClusters = "list",
         fun = "function"
         )
         )

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
## @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="compareClusterResult"),
          function (object){
              geneClusterLen <- length(object@geneClusters)
              cat ("Result of Comparing", geneClusterLen, "gene clusters", "\n")
          }
          )

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
## @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="compareClusterResult"),
          function(object) {
              return(object@compareClusterResult)
          }
          )

##' @rdname plot-methods
##' @aliases plot,compareClusterResult,ANY-method
##' @param x compareClusterResult object
##' @param type one of bar or dot
##' @param colorBy one of pvalue or p.adjust
##' @param showCategory category numbers
##' @param by one of geneRatio, Percentage or count
##' @param includeAll logical 
##' @param font.size font size
##' @param title figure title
setMethod("plot", signature(x="compareClusterResult"),
          function(x,
                   type="dot",
                   colorBy="p.adjust",
                   showCategory=5,
                   by="geneRatio",
                   includeAll=TRUE,
                   font.size=12,
                   title=""
                   ) {
              if (type == "dot" || type == "dotplot") {
                  dotplot(x, colorBy, showCategory, by, includeAll, font.size, title)
              } else if (type == "bar" || type == "barplot") {
                  barplot.compareClusterResult(x, colorBy, showCategory, by, includeAll, font.size, title)
              } else {
                  stop("type should be one of 'dot' or 'bar'...")
              }              
          })
##' dot plot method
##'
##'
##' @docType methods
##' @title dotplot
##' @rdname dotplot-methods
##' @aliases dotplot,compareClusterResult,ANY-method
##' @param object compareClusterResult object
##' @param colorBy one of pvalue or p.adjust
##' @param showCategory category numbers
##' @param by one of geneRatio, Percentage or count
##' @param includeAll logical 
##' @param font.size font size
##' @param title figure title
##' @importFrom DOSE dotplot
##' @exportMethod dotplot
setMethod("dotplot", signature(object="compareClusterResult"),
          function(object,
                   colorBy="p.adjust",
                   showCategory=5,
                   by="geneRatio",
                   includeAll=TRUE,
                   font.size=12,
                   title=""
                   ) {
              dotplot.compareClusterResult(object, colorBy, showCategory, by, includeAll, font.size, title)
          })


barplot.compareClusterResult <- function(height, colorBy="p.adjust", showCategory=5,
                                         by="geneRatio", includeAll=TRUE, font.size=12, title="", ...) {
    ## use *height* to satisy barplot generic definition
    ## actually here is an compareClusterResult object.
    df <- fortify(height, showCategory=showCategory, by=by, includeAll=includeAll)
    plotting.clusterProfile(df, type="bar", colorBy=colorBy, by=by, title=title, font.size=font.size)
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
    x <- lapply(enrichResultList, summary)
    names(x) <- names(enrichResultList)
    y <- ldply(x, "rbind")
    y <- rename(y, c(.id="Cluster"))
    y$Cluster = factor(y$Cluster, levels=names(enrichResultList))
    new("compareClusterResult",
        compareClusterResult = y)
    
}
