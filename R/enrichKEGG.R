enrichKEGG <- function(gene, organism="human", pvalueCutoff = 0.01) {
	pathID2ExtID <- as.list(KEGGPATHID2EXTID)
	#pathID <- mappedkeys(KEGGPATHID2EXTID)
	pathID <- names(pathID2ExtID)
	if (organism == "human") {
		idx <- grep("^hsa", pathID) ## select human pathways.
	} else if (organism == "mouse") {
		idx <- grep("^mmu", pathID) ## select mouse pathways.
	} else if (organism == "yeast") {
		idx <- grep("^sce", pathID) ## select yeast pathways.
	} else {
		stop (" Not supported yet... \n" )
	}
	orgPath2ExtID <- pathID2ExtID[idx]
	orgExtID <- unique(unlist(orgPath2ExtID))
		
	M = sapply(orgPath2ExtID, length)
	geneID.list = lapply(orgPath2ExtID, function(i) gene[gene %in% i])
	geneID <- sapply(geneID.list, function(i) paste(i, collapse="/"))
	k = sapply(geneID.list, length)
	#k = sapply(orgPath2ExtID, function(i) sum(gene %in% i))
	pathNum <- length(M)
	N <- rep(length(orgExtID), pathNum)
	n <- rep(length(gene), pathNum)
	args.df <- data.frame(numWdrawn=k-1, numW=M, numB=N-M, numDrawn=n)

	pvalues <- mdply(args.df, HyperG)
	pvalues <- pvalues[,5]
	
	GeneRatio <- mdply(data.frame(a=k, b=n), .yPaste)
	GeneRatio <- GeneRatio[,3]
	BgRatio <- mdply(data.frame(a=M, b=N), .yPaste)
	BgRatio <- BgRatio[,3]	
	pathwayID <- names(orgPath2ExtID)
	Description <- unlist(path2Name(pathwayID))
	
	keggOver <- data.frame(pathwayID=pathwayID, Description=Description, GeneRatio=GeneRatio, BgRatio=BgRatio, pvalue=pvalues)
	
	#qvalue =  fdrtool(keggOver$pvalue, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
	qobj = qvalue(keggOver$pvalue)
	qvalues <- qobj$qvalues
	keggOver <- data.frame(keggOver, qvalue=qvalues, geneID=geneID, Count=k)
	keggOver <- keggOver[order(pvalues),]
	keggOver <- keggOver[ keggOver$pvalue <= pvalueCutoff, ]
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
	function(x, caption="", font.size=12) {
		enrichKEGGResult <- x@enrichKEGGResult
		p <- .barplotInternal(enrichKEGGResult, caption, font.size)
		###color scale based on pvalue
		p + aes(fill=pvalue)
	}
)

