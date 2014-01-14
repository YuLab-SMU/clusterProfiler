##' Compare gene clusters functional profile
##' Given a list of gene set, this function will compute profiles of each gene
##' cluster.
##'
##'
##' @param geneClusters a list of entrez gene id.
##' @param fun One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .
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
##' 	data(gcSample)
##' 	xx <- compareCluster(gcSample, fun="enrichKEGG", organism="human", pvalueCutoff=0.05)
##' 	#summary(xx)
##' 	#plot(xx, type="dot", caption="KEGG Enrichment Comparison")
##'
compareCluster <- function(geneClusters, fun="enrichGO", ...) {
    fun <- eval(parse(text=fun))
    clProf <- llply(geneClusters,
                    .fun=function(i) {
			x=fun(i, ...)
                        if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
                            summary(x)
                        }
                    }
                    )

    clProf.df <- ldply(clProf, rbind)
    ##colnames(clProf.df)[1] <- "Cluster"
    clProf.df <- rename(clProf.df, c(.id="Cluster"))
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

##' show method for \code{compareClusterResult} instance
##'
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{compareClusterResult} instance.
##' @return message
##' @importFrom methods show
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="compareClusterResult"),
          function (object){
              geneClusterLen <- length(object@geneClusters)
                                        #fun <- object@fun
                                        #fun <- as.character(substitute(fun))
                                        #if (fun == "enrichKEGG") {
                                        #                analysis <- "KEGG Enrichment Analysis"
                                        #            } else if (fun == "groupGO") {
                                        #                analysis <- "GO Profiling Analysis"
                                        #           } else if (fun == "enrichGO") {
                                        #                analysis <- "GO Enrichment Analysis"
                                        #           } else if (fun == "enrichDO") {
                                        #              analysis <- "DO Enrichment Analysis"
                                        #	      } else {
                                        #		analysis <- "User specify Analysis"
                                        #	      }
                                        #              cat ("Compare", geneClusterLen, "gene clusters using", analysis, "\n")
              cat ("Result of Comparing", geneClusterLen, "gene clusters", "\n")
          }
          )

##' summary method for \code{compareClusterResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{compareClusterResult} instance.
##' @return A data frame
##' @importFrom stats4 summary
##' @exportMethod summary
##' @author Guangchuang Yu \url{http://ygc.name}
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
setMethod("plot", signature(x="compareClusterResult"),
          function(x,
                   type="dot",
                   title="",
                   font.size=12,
                   showCategory=5,
                   by="geneRatio",
                   colorBy="p.adjust",
                   angle.axis.x=90,
                   hjust.axis.x=1,
                   vjust.axis.x=0.5) {

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
                  ## first try
                  ## cls <- as.character(result$Cluster)
                  ## clsu <- unique(cls)
                  ## idx <- sapply(clsu, function(i) which(i == cls)[1])
                  ## gcSize <- result$Count[idx]
                  ## names(gcSize) <- clsu
                  ## result$GeneRatio <- result$Count / gcSize[cls]
                  ## result$Cluster <- paste(cls, "(", gcSize[cls], ")", sep="")

                  ## second try
                  ## result$GeneRatio <- sapply(as.character(result$GeneRatio), function(i) eval(parse(text=i)))

                  ## final way
                  gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
                  gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
                  result$GeneRatio = gsize/gcsize
              } else {
                  ## nothing
              }
              p <- plotting.clusterProfile(result,
                                           type=type,
                                           by=by,
                                           colorBy=colorBy,
                                           title=title,
                                           font.size=font.size,
                                           angle.axis.x=angle.axis.x,
                                           hjust.axis.x=hjust.axis.x,
                                           vjust.axis.x=vjust.axis.x)
              return(p)
          }
          )
