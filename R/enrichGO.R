

#' GO Enrichment Analysis of a gene set.
#' Given a vector of genes, this function will return the enrichment GO
#' categories with FDR control.
#' 
#' 
#' @param gene a vector of entrez gene id.
#' @param organism Currently, only "human", "mouse" and "yeast" supported.
#' @param ont One of "MF", "BP", and "CC" subontologies.
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param readable if readable is TRUE, the gene IDs will mapping to gene
#'   symbols.
#' @return A \code{enrichGOResult} instance.
#' @seealso \code{\link{enrichGOResult-class}}, \code{\link{compareCluster}}
#' @keywords manip
#' @examples
#' 
#' 	#data(gcSample)
#' 	#yy <- enrichGO(gcSample[[1]], organism="human", ont="BP", pvalueCutoff=0.01, testDirection="over")
#' 	#head(summary(yy))
#' 	#plot(yy)
#' 
enrichGO <- function(gene, organism="human", ont="MF", pvalueCutoff=0.01, readable=FALSE) {
    goterms <- Ontology(GOTERM)
    goterms <- names(goterms[goterms == ont])

    orgTerm <- switch(organism,
                      human = mappedkeys(org.Hs.egGO2ALLEGS),
                      mouse = mappedkeys(org.Mm.egGO2ALLEGS),
                      yeast = mappedkeys(org.Sc.sgdGO2ALLORFS),
                      )

    Terms <- goterms[goterms %in% orgTerm]

    GO2ExtID <- switch(organism,
                       human = mget(Terms, org.Hs.egGO2ALLEGS, ifnotfound=NA),
                       mouse = mget(Terms, org.Mm.egGO2ALLEGS, ifnotfound=NA),
                       yeast = mget(Terms, org.Sc.sgdGO2ALLORFS, ifnotfound=NA),
                       )
	GO2ExtID <- lapply(GO2ExtID, function(i) unique(i))

    orgExtID <- unique(unlist(GO2ExtID))

    geneID.list = lapply(GO2ExtID, function(i) gene[gene %in% i])
    if (readable) {
        geneID.list <- geneID2geneName(geneID.list, organism)
    }
    geneID <- sapply(geneID.list, function(i) paste(i, collapse="/"))


    k = sapply(geneID.list, length)
    retain <- which(k != 0)
	k = k[retain]
	GO2ExtID <- GO2ExtID[retain]
	geneID=geneID[retain]

	M <- sapply(GO2ExtID, length)

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
    GOID <- names(GO2ExtID)
    Description <- unlist(sapply(GOID, .GO2Term))

    goOver <- data.frame(GOID=GOID, Description=Description, GeneRatio=GeneRatio, BgRatio=BgRatio, pvalue=pvalues)

                                        #qvalue =  fdrtool(keggOver$pvalue, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
    qobj = qvalue(goOver$pvalue, lambda=0.05, pi0.method="bootstrap")
    qvalues <- qobj$qvalues
    goOver <- data.frame(goOver, qvalue=qvalues, geneID=geneID, Count=k)
    goOver <- goOver[order(goOver$pvalue),]
    goOver <- goOver[ goOver$pvalue <= pvalueCutoff, ]
    goOver$Description <- as.character(goOver$Description)

    new("enrichGOResult",
        enrichGOResult = goOver,
        pvalueCutoff=pvalueCutoff,
        Organism = organism,
        Ont = ont,
        Gene = gene
        )
}




#enrichGO <- function(gene, organism="human", ont="MF", pvalueCutoff = 0.01, testDirection="over") {
#	if (organism == "human") {
#		annotation="org.Hs.eg.db"
#		geneUniverse =  mappedkeys(org.Hs.egGO)
#	} else if (organism == "mouse") {
#		annotation="org.Mm.eg.db"
#		geneUniverse =  mappedkeys(org.Mm.egGO)
#	} else if (organism == "yeast"){
#		annotation="org.Sc.sgd.db"
#		geneUniverse =  mappedkeys(org.Sc.sgdGO)
#	}else {
#		stop (" Not supported yet... \n" )
#	}
#	params = new("GOHyperGParams", geneIds=gene, universeGeneIds=geneUniverse, annotation=annotation, ontology=ont, pvalueCutoff = 1, conditional = FALSE, testDirection=testDirection)
#	hgOver = hyperGTest(params)
#	hgOver.df <- summary(hgOver)
#	if (nrow(hgOver.df) == 0) {
#		return(NA)
#	}
#	colnames(hgOver.df)[c(1,7)] <- c("GOID", "Description")
#
#	## get GeneIDs annotated with GOID.
#	goGene <- .goGene(GOID=hgOver.df[,1], gene, organism)
#	GeneIDs = .getGeneID(goGene)
#
#	GeneSetSize <- rep(length(hgOver@geneIds), nrow(hgOver.df))
#
#	#qvalue =  fdrtool(hgOver.df$Pvalue, statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
#	#p.adjust(hgOver.df$Pvalue, method='fdr')
#	qobj <- qvalue(hgOver.df$Pvalue)
#	qvalues <- qobj$qvalues
#	hgOver.df <- data.frame(hgOver.df, GeneSetSize=GeneSetSize, GeneID=GeneIDs, qvalue=qvalues)
#	hgOver.df <- hgOver.df[,c(1,7,2,10,3:5,8,6,9)]
#
#	hgOver.df <- hgOver.df[ hgOver.df$Pvalue <= pvalueCutoff, ]
#
#	new("enrichGOResult",
#		enrichGOResult = hgOver.df,
#		pvalueCutoff=pvalueCutoff,
#		testDirection=testDirection,
#		Ont = ont,
#		Organism = organism,
#		Gene = gene
#	)
#}

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
          function(x, caption="", font.size=12) {
              enrichGOResult <- summary(x)
              p <- .barplotInternal(enrichGOResult, caption, font.size)
###color scale based on pvalue
              p + aes(fill=Pvalue)
          }
          )

#setMethod("plot", signature(x="GOHyperGResult"),
#	function(x, caption="", font.size=12) {
#		enrichGOResult <- summary(x)
#		#colnames(enrichGOResult)[7] <- "Description"
#		enrichGOResult <- rename(enrichGOResult, c(Term="Description"))
#		p <- .barplotInternal(enrichGOResult, caption, font.size)
#		###color scale based on pvalue
#		p + aes(fill=Pvalue)
#	}
#)
