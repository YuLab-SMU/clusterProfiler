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
##' 	data(gcSample)
##' 	xx <- compareCluster(gcSample, fun="enrichKEGG", organism="human", pvalueCutoff=0.05)
##' 	#summary(xx)
##' 	#plot(xx, type="dot", caption="KEGG Enrichment Comparison")
##'
compareCluster <- function(geneClusters, fun="enrichGO", data='', ...) {

    fun <- eval(parse(text=fun))
    # Use formula interface for compareCluster
    if (typeof(geneClusters) == 'language') {
        if (!is.data.frame(data)) {
            print ('no data provided with formula')
        } else {
#            data[[all.vars(geneClusters)[1]]] = as.character(data[[all.vars(geneClusters)[1]]])
##            geneClusters = by(data[all.vars(geneClusters)[1]], data[all.vars(geneClusters)[2]], function(x) list(as.character(x)))

#            geneClusters = split(data[all.vars(geneClusters)[1]], data[all.vars(geneClusters)[2]])
            geneClusters = dlply(.data=data, all.vars(geneClusters)[2], .fun=function(x) {as.character(x[all.vars(geneClusters)[1]])})
#            clProf <- dlply(geneClusters, .fun=function(i))
#            clProf.df <- aggregate(geneClusters, data, function(i) { print(head(i)); x=fun(as.character(i))
#                            if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
#                                summary(x)
#                            }
#                        })
#            geneClusters = unique(data[all.vars(geneClusters)[1]])
        }   
        print('dsadasd')
#        print(clProf.df)
#        print(summary(clProf.df))
    }
#    } else {

        print(summary(geneClusters))
        print('dasda')
        clProf <- llply(geneClusters,
                        .fun=function(i) {
                x=fun(i, ...)
                            if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
                                summary(x)
                            }
                        }
                        )
#    }
    print(summary(clProf))
    clProf.df <- ldply(clProf, rbind)
    print(head(clProf.df %>% dplyr::select(-geneID)))
    clProf.df <- rename(clProf.df, c(.id="Cluster"))

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
setMethod("plot", signature(x="compareClusterResult"),
          function(x,
                   type="dot",
                   title="",
                   font.size=12,
                   showCategory=5,
                   by="geneRatio",
                   colorBy="p.adjust") {

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
                  result$Cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")   # Transform to factor to keep the original order?
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
