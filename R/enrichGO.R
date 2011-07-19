

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
#' @export
#' @examples
#'
#' 	#data(gcSample)
#' 	#yy <- enrichGO(gcSample[[1]], organism="human", ont="BP", pvalueCutoff=0.01)
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

    GO2ExtID <- getGO2ExtID(Terms, organism)

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


    GeneRatio <- mdply(data.frame(a=k, b=n), getRatio)
    GeneRatio <- GeneRatio[,3]
    BgRatio <- mdply(data.frame(a=M, b=N), getRatio)
    BgRatio <- BgRatio[,3]
    GOID <- names(GO2ExtID)
    Description <- GO2Term(GOID)

    goOver <- data.frame(GOID=GOID, Description=Description, GeneRatio=GeneRatio, BgRatio=BgRatio, pvalue=pvalues)

    qobj = qvalue(goOver$pvalue, lambda=0.05, pi0.method="bootstrap")
    qvalues <- qobj$qvalues
    goOver <- data.frame(goOver, qvalue=qvalues, geneID=geneID, Count=k)
    goOver <- goOver[order(goOver$pvalue),]
    goOver <- goOver[ goOver$pvalue <= pvalueCutoff, ]
    goOver$Description <- as.character(goOver$Description)
    rownames(goOver) <- goOver$GOID

    new("enrichGOResult",
        enrichGOResult = goOver,
        pvalueCutoff=pvalueCutoff,
        Organism = organism,
        Ont = ont,
        Gene = gene
        )
}

##' An S4 class that stores Gene Ontology enrichment result
##' @slot enrichGOResult GO enrichment result
##' @slot pvalueCutoff pvalueCutoff
##' @slot Ont Ontology
##' @slot Organism one of "human", "mouse" and "yeast"
##' @slot Gene Gene IDs
##' @author Guangchuang Yu
setClass("enrichGOResult",
         representation=representation(
         enrichGOResult="data.frame",
         pvalueCutoff="numeric",
         Ont = "character",
         Organism = "character",
         Gene = "character"
         )
         )

##' show method for \code{enrichGOResult} instance
##'
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{enrichGOResult} instance.
##' @return message
##' @author GuangchuangYu
setMethod("show", signature(object="enrichGOResult"),
          function (object){
              ont = object@Ont
              Organism = object@Organism
              GeneNum = length(object@Gene)
              pvalueCutoff=object@pvalueCutoff
              cat ("Hypergeometric test of over-representation GO (", ont, ") categories", " for ", GeneNum, Organism, "genes\n", "p value <", pvalueCutoff, "\n")
          }
          )

##' summary method for \code{enrichGOResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{enrichGOResult} instance.
##' @return A data frame
##' @author GuangchuangYu
setMethod("summary", signature(object="enrichGOResult"),
          function(object) {
              return(object@enrichGOResult)
          }
          )

##' plot method for \code{enrichGOResult} instance
##'
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##'
##' @title plot method
##' @param x A \code{enrichGOResult} instance.
##' @param title graph title
##' @param font.size graph font size
##' @return ggplot object
##' @author Guangchuang Yu
setMethod("plot", signature(x="enrichGOResult"),
          function(x, title="", font.size=12) {
              enrichGOResult <- summary(x)
              p <- .barplotInternal(enrichGOResult, title, font.size)
              ##color scale based on pvalue
              p + aes(fill=Pvalue)
          }
          )
