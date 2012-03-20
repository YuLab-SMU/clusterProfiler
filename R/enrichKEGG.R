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
##' @return A \code{enrichResult} instance.
##' @export
##' @importFrom DOSE enrich.internal
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @importMethodsFrom DOSE setReadable
##' @importFrom DOSE EXTID2NAME


##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @importClassesFrom methods data.frame
##' @importFrom KEGG.db KEGGPATHID2EXTID


##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @examples
##'
##' 	data(gcSample)
##' 	yy = enrichKEGG(gcSample[[5]], pvalueCutoff=0.01)
##' 	head(summary(yy))
##' 	#plot(yy)
##'
enrichKEGG <- function(gene,
                       organism="human",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable=FALSE) {

    enrich.internal(gene,
                    organism = organism,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    ont = "KEGG",
                    readable = readable)

}


##' @importFrom DOSE EXTID2TERMID
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom KEGG.db KEGGEXTID2PATHID
##' @S3method EXTID2TERMID KEGG
EXTID2TERMID.KEGG <- function(gene, organism) {
    gene <- as.character(gene)
    qExtID2PathID <- mget(gene, KEGGEXTID2PATHID, ifnotfound=NA)
    notNA.idx <- unlist(lapply(qExtID2PathID, function(i) !all(is.na(i))))
    qExtID2PathID <- qExtID2PathID[notNA.idx]
    return(qExtID2PathID)
}

##' @importFrom DOSE TERMID2EXTID
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom KEGG.db KEGGPATHID2EXTID
##' @S3method TERMID2EXTID KEGG
TERMID2EXTID.KEGG <- function(term, organism) {
    pathID2ExtID <- mget(unique(term), KEGGPATHID2EXTID, ifnotfound=NA)
    return(pathID2ExtID)
}

##' @importFrom DOSE ALLEXTID
##' @importFrom KEGG.db KEGGPATHID2EXTID
##' @S3method ALLEXTID KEGG
ALLEXTID.KEGG <- function(organism) {
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
    class(orgPathID) <- "KEGG"
    orgPath2ExtID <- TERMID2EXTID(orgPathID)
    orgPath2ExtID <- lapply(orgPath2ExtID, function(i) unique(i))

    orgExtID <- unique(unlist(orgPath2ExtID))
    return(orgExtID)
}

##' @importFrom DOSE TERM2NAME
##' @importFrom KEGG.db KEGGPATHID2NAME
##' @importMethodsFrom AnnotationDbi mget
##' @S3method TERM2NAME KEGG
TERM2NAME.KEGG <- function(term) {
    term <- as.character(term)
    pathIDs <- gsub("^\\D+", "",term, perl=T)
    path2name <- unlist(mget(pathIDs, KEGGPATHID2NAME))
    return(path2name)
}
