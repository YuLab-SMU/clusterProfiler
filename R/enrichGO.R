enrichGO <- function(gene, organism="human", ont="MF", pvalueCutoff = 0.01, testDirection="over") {
	if (organism == "human") {
		#require(org.Hs.eg.db)
		annotation="org.Hs.eg.db"
		geneUniverse =  mappedkeys(org.Hs.egGO)
	} else if (organism == "mouse") {
		#require(org.Mm.eg.db)
		annotation="org.Mm.eg.db"
		geneUniverse =  mappedkeys(org.Mm.egGO)		
	} else {
		stop (" Not supported yet... \n" )
	}
	params = new("GOHyperGParams", geneIds=gene, universeGeneIds=geneUniverse, annotation=annotation, ontology=ont, pvalueCutoff = pvalueCutoff, conditional = FALSE, testDirection=testDirection)
	hgOver = hyperGTest(params)
	hgOver.df <- summary(hgOver)
	if (nrow(hgOver.df) == 0) {
		return(NA)
	}
	colnames(hgOver.df)[c(1,7)] <- c("GOID", "Description")
	
	## get GeneIDs annotated with GOID.
	goGene <- .goGene(GOID=hgOver.df[,1], gene, organism) 
	GeneIDs = .getGeneID(goGene)
	
	GeneSetSize <- rep(length(hgOver@geneIds), nrow(hgOver.df))
	
	qvalue =  fdrtool(hgOver.df$Pvalue, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval 	  
	#p.adjust(hgOver.df$Pvalue, method='fdr')

	hgOver.df <- data.frame(hgOver.df, GeneSetSize=GeneSetSize, GeneID=GeneIDs, qvalue=qvalue)
	hgOver.df <- hgOver.df[,c(1,7,2,10,3:5,8,6,9)]
	
	hgOver.df <- hgOver.df[ hgOver.df$Pvalue <= pvalueCutoff, ]
	
	new("enrichGOResult", 
		enrichGOResult = hgOver.df,
		pvalueCutoff=pvalueCutoff,
		testDirection=testDirection,
		Ont = ont,
		Organism = organism,
		Gene = gene
	)
}

setClass("enrichGOResult",
	representation=representation(
		enrichGOResult="data.frame",
		pvalueCutoff="numeric",
		testDirection="character",
		Ont = "character",
		Organism = "character",
		Gene = "character"
	)
)

setMethod("show", signature(object="enrichGOResult"),
	function (object){
		ont = object@Ont
		Organism = object@Organism
		GeneNum = length(object@Gene)
		testDirection = object@testDirection
		pvalueCutoff=object@pvalueCutoff
		cat (GeneNum, Organism, "Genes to GO", ont, "test for", paste(testDirection, "-representation.", sep=""), "\n", "p value <", pvalueCutoff, "\n")
	}
)

setMethod("summary", signature(object="enrichGOResult"),
	function(object) {
		return(object@enrichGOResult)
	}
)

setMethod("plot", signature(x="enrichGOResult"),
	function(x, caption="") {
		enrichGOResult <- summary(x)
		p <- .barplotInternal(enrichGOResult, caption)
		###color scale based on pvalue
		p + aes(fill=Pvalue)
	}
)

setMethod("plot", signature(x="GOHyperGResult"),
	function(x, caption="") {
		enrichGOResult <- summary(x)
		colnames(enrichGOResult)[7] <- "Description"
		p <- .barplotInternal(enrichGOResult, caption)
		###color scale based on pvalue
		p + aes(fill=Pvalue)
	}
)