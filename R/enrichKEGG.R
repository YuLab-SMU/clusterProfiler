##' KEGG Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment KEGG
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param organism Currently, only "human" and "mouse" supported.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param qvalueCutoff Cutoff value of qvalue.
##' @param readable if readable is TRUE, the gene IDs will mapping to gene
##'   symbols.
##' @return A \code{enrichKEGGResult} instance.
##' @export
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @importClassesFrom methods data.frame
##' @importFrom methods new
##' @importFrom qvalue qvalue
##' @importFrom KEGG.db KEGGPATHID2EXTID
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{enrichKEGGResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @examples
##'
##' 	data(gcSample)
##' 	yy = enrichKEGG(gcSample[[5]], pvalueCutoff=0.01)
##' 	head(summary(yy))
##' 	#plot(yy)
##'
enrichKEGG <- function(gene, organism="human",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                       readable=FALSE) {
    ##pathID2ExtID <- as.list(KEGGPATHID2EXTID)
    ##pathID <- names(pathID2ExtID)
    pathID <- mappedkeys(KEGGPATHID2EXTID)

    ## select species specific pathways
    if (organism == "human") {
        idx <- grep("^hsa", pathID)
    } else if (organism == "mouse") {
        idx <- grep("^mmu", pathID)
    } else if (organism == "yeast") {
        idx <- grep("^sce", pathID)
    } else {
        stop (" Not supported yet... \n" )
    }
    ##orgPath2ExtID <- pathID2ExtID[idx]
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
    ##pvalues <- mdply(args.df, HyperG)
    ##pvalues <- pvalues[,5]
    pvalues <- apply(args.df, 1, HyperG)

    ##GeneRatio <- mdply(data.frame(a=k, b=n), getRatio)
    ##GeneRatio <- GeneRatio[,3]
    GeneRatio <- apply(data.frame(a=k, b=n), 1, getRatio)

    ##BgRatio <- mdply(data.frame(a=M, b=N), getRatio)
    ##BgRatio <- BgRatio[,3]
    BgRatio <- apply(data.frame(a=M, b=N), 1, getRatio)

    pathwayID <- names(orgPath2ExtID)
    Description <- unlist(path2Name(pathwayID))

    keggOver <- data.frame(pathwayID=pathwayID,
                           Description=Description,
                           GeneRatio=GeneRatio,
                           BgRatio=BgRatio,
                           pvalue=pvalues)


    qobj = qvalue(keggOver$pvalue, lambda=0.05, pi0.method="bootstrap")
    qvalues <- qobj$qvalues
    keggOver <- data.frame(keggOver, qvalue=qvalues,
                           geneID=geneID, Count=k)
    keggOver <- keggOver[order(pvalues),]

    keggOver <- keggOver[ keggOver$pvalue <= pvalueCutoff, ]
    keggOver <- keggOver[ keggOver$qvalue <= qvalueCutoff, ]

    keggOver$Description <- as.character(keggOver$Description)

    new("enrichKEGGResult",
        enrichKEGGResult = keggOver,
        pvalueCutoff=pvalueCutoff,
        Organism = organism,
        Gene = gene
	)
}

##' Class "enrichKEGGResult"
##' This class represents the result of KEGG enrichment analysis.
##'
##'
##' @name enrichKEGGResult-class
##' @aliases enrichKEGGResult-class show,enrichKEGGResult-method
##'   summary,enrichKEGGResult-method plot,enrichKEGGResult-method
##' @docType class
##' @slot enrichKEGGResult KEGG enrichment result
##' @slot pvalueCutoff pvalueCutoff
##' @slot qvalueCutoff qvalueCutoff
##' @slot Organism one of "humna", "mouse", and "yeast"
##' @slot Gene Gene IDs
##' @exportClass enrichKEGGResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\linkS4class{compareClusterResult}}
##'   \code{\link{compareCluster}} \code{\link{enrichKEGG}}
##' @keywords classes
setClass("enrichKEGGResult",
         representation=representation(
         enrichKEGGResult="data.frame",
         pvalueCutoff="numeric",
         qvalueCutoff="numeric",
         Organism = "character",
         Gene = "character"
         )
         )

##' show method for \code{enrichKEGGResult} instance
##'
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{enrichKEGGResult} instance.
##' @return message
##' @importFrom methods show
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="enrichKEGGResult"),
          function (object){
              Organism = object@Organism
              GeneNum = length(object@Gene)
              pvalueCutoff=object@pvalueCutoff
              cat (GeneNum, Organism,
                   "Genes to KEGG test for over-representation.", "\n",
                   "p value <", pvalueCutoff, "\n")
          }
          )

##' summary method for \code{enrichKEGGResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{enrichKEGGResult} instance.
##' @return A data frame
##' @importFrom stats4 summary
##' @exportMethod summary
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="enrichKEGGResult"),
          function(object) {
              return(object@enrichKEGGResult)
          }
          )

##' plot method for \code{enrichKEGGResult} instance
##'
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##'
##' @title plot method
##' @param x A \code{enrichKEGGResult} instance.
##' @param title graph title
##' @param font.size graph font size
##' @param showCategory number of KEGG categories to show.
##' @return ggplot object
##' @importFrom graphics plot
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 scale_fill_continuous
##' @exportMethod plot
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="enrichKEGGResult"),
          function(x, title="", font.size=12, showCategory=5) {
              enrichKEGGResult <- x@enrichKEGGResult
              if ( is.numeric(showCategory) & showCategory < nrow(enrichKEGGResult) ) {
                  enrichKEGGResult <- enrichKEGGResult[1:showCategory,]
              }
              p <- plotting.barplot(enrichKEGGResult, title, font.size)
              ##color scale based on pvalue
              p <- p +
                  aes(fill=pvalue) +
                      scale_fill_continuous(low="red", high="blue")
              return(p)
          }
          )

