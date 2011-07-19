

#' Compare gene clusters functional profile
#' Given a list of gene set, this function will compute profiles of each gene
#' cluster.
#'
#'
#' @param geneClusters a list of entrez gene id.
#' @param fun One of groupGO and enrichGO.
#' @param ...  Other arguments.
#' @return A \code{clusterProfResult} instance.
#' @seealso \code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
#'   \code{\link{enrichGO}}
#' @keywords manip
#' @examples
#'
#' 	data(gcSample)
#' 	xx <- compareCluster(gcSample, fun=enrichKEGG, organism="human", pvalueCutoff=0.05)
#' 	#summary(xx)
#' 	#plot(xx, type="dot", caption="KEGG Enrichment Comparison")
#'


#' Compare gene clusters functional profile
#' Given a list of gene set, this function will compute profiles of each gene
#' cluster.
#' 
#' 
#' @param geneClusters a list of entrez gene id.
#' @param fun One of groupGO and enrichGO.
#' @param ...  Other arguments.
#' @return A \code{clusterProfResult} instance.
#' @seealso \code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
#'   \code{\link{enrichGO}}
#' @keywords manip
#' @examples
#' 
#' 	data(gcSample)
#' 	xx <- compareCluster(gcSample, fun=enrichKEGG, organism="human", pvalueCutoff=0.05)
#' 	#summary(xx)
#' 	#plot(xx, type="dot", caption="KEGG Enrichment Comparison")
#' 
compareCluster <- function(geneClusters, fun=enrichGO, ...) {
    clProf <- llply(geneClusters,
                    .fun=function(i) {
			x=fun(i, ...)
			summary(x)
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

setClass("compareClusterResult",
         representation = representation(
         compareClusterResult = "data.frame",
         geneClusters = "list",
         fun = "function"
         )
         )

setMethod("show", signature(object="compareClusterResult"),
          function (object){
              geneClusterLen <- length(object@geneClusters)
              fun <- object@fun
              fun <- as.character(substitute(fun))
              if (fun == "enrichKEGG") {
                  analysis <- "KEGG Enrichment Analysis"
              } else if (fun == "groupGO") {
                  analysis <- "GO Profiling Analysis"
              } else if (fun == "enrichGO") {
                  analysis <- "GO Enrichment Analysis"
              }
              cat ("Compare", geneClusterLen, "gene clusters using", analysis, "\n")
          }
          )

setMethod("summary", signature(object="compareClusterResult"),
          function(object) {
              return(object@compareClusterResult)
          }
          )

setMethod("plot", signature(x="compareClusterResult"),
          function(x, type="dot", title="", font.size=12, limit=5, by="percentage") {
              clProf.df <- summary(x)

              ## get top 5 (default) categories of each gene cluster.
              if (is.null(limit)) {
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
                                  N=limit
                                  )
              }

              ## remove zero count
              result$Description <- as.character(result$Description) ## un-factor
              GOlevel <- result[,c(2,3)] ## GO ID and Term
              GOlevel <- unique(GOlevel)

              result <- result[result$Count != 0, ]
              result$Description <- factor(result$Description, levels=rev(GOlevel[,2]))


              if (by=="percentage") {
                  Description <- Count <- NULL # to satisfy codetools
                  result <- ddply(result, .(Description), transform, Percentage = Count/sum(Count), Total = sum(Count))

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
                  ##result[, colnames(result) != "Total"]

                  result$Description <- factor(result$Description, levels=rev(Termlevel))

              } else if (by == "count") {

              } else {

              }
              .PlotClusterProfInternal(result, type, by, title, font.size)
          }
          )
