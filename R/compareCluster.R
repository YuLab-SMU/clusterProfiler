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
##' @importFrom plyr ddply
##' @importFrom plyr mdply
##' @importFrom plyr .
##' @param x compareClusterResult object
##' @param type one of bar or dot
##' @param title figure title
##' @param font.size font size
##' @param showCategory category numbers
##' @param by one of geneRatio, Percentage or count
##' @param colorBy one of pvalue or p.adjust
##' @param includeAll logical 
setMethod("plot", signature(x="compareClusterResult"),
          function(x,
                   type="dot",
                   title="",
                   font.size=12,
                   showCategory=5,
                   by="geneRatio",
                   colorBy="p.adjust",
                   includeAll=TRUE
                   ) {

              clProf.df <- summary(x)

              ## get top 5 (default) categories of each gene cluster.
              if (is.null(showCategory)) {
                  result <- clProf.df
              } else {
                  Cluster <- NULL # to satisfy codetools
                  result <- ddply(.data = clProf.df,
                                  .variables = .(Cluster),
                                  .fun = function(df, N) {
                                      if (length(df$Count) > N) {
                                          idx <- order(df$Count, decreasing=T)[1:N]
                                          return(df[idx,])
                                      } else {
                                          return(df)
                                      }
                                  },
                                  N=showCategory
                                  )

              }
              if (includeAll == TRUE) {
                  result = subset(clProf.df, ID %in% result$ID)
              }

              ## remove zero count
              result$Description <- as.character(result$Description) ## un-factor
              GOlevel <- result[,c(2,3)] ## GO ID and Term
              GOlevel <- unique(GOlevel)

              result <- result[result$Count != 0, ]
              result$Description <- factor(result$Description,
                                           levels=rev(GOlevel[,2]))


              if (by=="rowPercentage") {
                  Description <- Count <- NULL # to satisfy codetools
                  result <- ddply(result,
                                  .(Description),
                                  transform,
                                  Percentage = Count/sum(Count),
                                  Total = sum(Count))

                  ## label GO Description with gene counts.
                  x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
                  y <- sapply(x[,3], paste, ")", sep="")
                  result$Description <- y

                  ## restore the original order of GO Description
                  xx <- result[,c(2,3)]
                  xx <- unique(xx)
                  rownames(xx) <- xx[,1]
                  Termlevel <- xx[as.character(GOlevel[,1]),2]

                  ##drop the *Total* column
                  result <- result[, colnames(result) != "Total"]

                  result$Description <- factor(result$Description,
                                               levels=rev(Termlevel))

              } else if (by == "count") {
                  ## nothing
              } else if (by == "geneRatio") {
                  gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
                  gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
                  result$GeneRatio = gsize/gcsize
                  result$Cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")
              } else {
                  ## nothing
              }
              p <- plotting.clusterProfile(result,
                                           type=type,
                                           by=by,
                                           colorBy=colorBy,
                                           title=title,
                                           font.size=font.size)
              return(p)
          }
          )
