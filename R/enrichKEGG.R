enrichKEGG <- function(gene, organism="human", pvalueCutoff = 0.05, readable=FALSE) {
                                        #pathID2ExtID <- as.list(KEGGPATHID2EXTID)
                                        #pathID <- names(pathID2ExtID)
    pathID <- mappedkeys(KEGGPATHID2EXTID)

    if (organism == "human") {
        idx <- grep("^hsa", pathID) ## select human pathways.
    } else if (organism == "mouse") {
        idx <- grep("^mmu", pathID) ## select mouse pathways.
    } else if (organism == "yeast") {
        idx <- grep("^sce", pathID) ## select yeast pathways.
    } else {
        stop (" Not supported yet... \n" )
    }
                                        #orgPath2ExtID <- pathID2ExtID[idx]
    orgPathID <- pathID[idx]
    orgPath2ExtID <- mget(orgPathID, KEGGPATHID2EXTID, ifnotfound=NA)
    orgPath2ExtID <- lapply(orgPath2ExtID, function(i) unique(i))

    orgExtID <- unique(unlist(orgPath2ExtID))

    geneID.list = lapply(orgPath2ExtID, function(i) gene[gene %in% i])
    if (readable) {
        geneID.list <- geneID2geneName(geneID.list, organism)
    }
    geneID <- sapply(geneID.list, function(i) paste(i, collapse="/"))
    k = sapply(geneID.list, length)
    retain <- which(k != 0)
    k = k[retain]
    orgPath2ExtID <- orgPath2ExtID[retain]
    geneID=geneID[retain]

    M = sapply(orgPath2ExtID, length)

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
    qobj = qvalue(keggOver$pvalue, lambda=0.05, pi0.method="bootstrap")
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

