enrichKEGG <- function(gene, organism="human", pvalueCutoff = 0.01) {
	if (organism == "human") {
		utils::data(list="hsa_ncbi-geneid", package="SubpathwayMiner")
	} else if (organism == "mouse") {
		utils::data(list="mmu_ncbi-geneid", package="clusterProfiler")
	} else {
		stop (" Not supported yet... \n" )
	}
	ann <- getAnn(unique(gene))
	geneIDs <- c()
	for (i in 1:length(ann)) {
		geneID <- ann[[i]][[2]]
		x = paste(geneID, collapse="/")
		geneIDs <- c(geneIDs, x)
	}
	keggOver <- printAnn(ann)
	Count <- sapply(keggOver$annGeneRatio, function(i) unlist(strsplit(as.character(i), split="/"))[1])
	Count <- as.numeric(Count)
	
	keggOver <- data.frame(pathwayID=rownames(keggOver), keggOver, geneID=geneIDs, Count=Count)
	colnames(keggOver)[2] <- "Description"
	idx = as.numeric(as.character(keggOver$pvalue)) < pvalueCutoff
	keggOver = keggOver[idx,]
	keggOver$pvalue <- as.numeric(as.character(keggOver$pvalue))
	keggOver$Description <- as.character(keggOver$Description)
	
	new("enrichKEGGResult", 
		enrichKEGGResult = keggOver,
		pvalueCutoff=pvalueCutoff,
		Organism = organism,
		Gene = gene
	)
}

setClass("enrichKEGGResult",
	representation=representation(
		enrichKEGGResult="data.frame",
		pvalueCutoff="numeric",
		Organism = "character",
		Gene = "character"
	)
)

setMethod("show", signature(object="enrichKEGGResult"),
	function (object){
		Organism = object@Organism
		GeneNum = length(object@Gene)
		pvalueCutoff=object@pvalueCutoff
		cat (GeneNum, Organism, "Genes to KEGG test for over-representation.", "\n", "p value <", pvalueCutoff, "\n")
	}
)

setMethod("summary", signature(object="enrichKEGGResult"),
	function(object) {
		return(object@enrichKEGGResult)
	}
)

setMethod("plot", signature(x="enrichKEGGResult"),
	function(x, caption="") {
		enrichKEGGResult <- x@enrichKEGGResult
		p <- .barplotInternal(enrichKEGGResult, caption)
		###color scale based on pvalue
		p + aes(fill=pvalue)
	}
)

