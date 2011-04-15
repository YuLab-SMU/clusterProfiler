compareCluster <- function(geneClusters, fun=enrichGO, ...) {
	clProf <- llply(geneClusters, 
		.fun=function(i) {
			x=fun(i, ...) 
			summary(x)
		}
	)
	
	clProf.df <- ldply(clProf, rbind)
	#colnames(clProf.df)[1] <- "Cluster"
	clProf.df <- rename(clProf.df, c(.id="Cluster"))
	new("compareClusterResult", 
		compareClusterResult = clProf.df,
		geneClusters = geneClusters
	)
}

setClass("compareClusterResult",
	representation=representation(
		compareClusterResult="data.frame",
		geneClusters = "list"
	)
)

setMethod("show", signature(object="compareClusterResult"),
	function (object){
		geneClusterLeng <- length(object@geneClusters)
		clProf.df <- object@compareClusterResult
		nc <- ncol(clProf.df)
		if (nc == 9) {
			analysis <- "KEGG Enrichment Analysis"
		} else if (nc == 5) {
			analysis <- "GO Profiling Analysis"
		} else {
			analysis <- "GO Enrichment Analysis"
		}
		cat ("Compare", geneClusterLeng, "gene clusters using", analysis, "\n")  
	}
)

setMethod("summary", signature(object="compareClusterResult"),
	function(object) {
		return(object@compareClusterResult)
	}
)

setMethod("plot", signature(x="compareClusterResult"),
	function(x, type="dot", caption="", limit=5, by="percentage") {
		clProf.df <- .compareClusterResultTopN(x, limit)
		clProf.df <- .removeZeroCount(clProf.df)	
		
		if (by=="percentage") {
			clProf.df <- .compareClusterResultReshape(clProf.df)
		} else if (by == "count") {
			
		} else {
		
		}
		.PlotClusterProfInternal(clProf.df, type, by, caption)
	}
)