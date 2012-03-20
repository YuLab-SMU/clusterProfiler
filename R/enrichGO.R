##' GO Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment GO
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param organism Currently, only "human", "mouse" and "yeast" supported.
##' @param ont One of "MF", "BP", and "CC" subontologies.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param qvalueCutoff Cutoff value of qvalue.
##' @param readable if readable is TRUE, the gene IDs will mapping to gene
##'   symbols.
##' @return A \code{enrichResult} instance.
##' @importFrom DOSE enrich.internal
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @importMethodsFrom DOSE setReadable
##' @importFrom DOSE EXTID2NAME
##' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @examples
##'
##' 	#data(gcSample)
##' 	#yy <- enrichGO(gcSample[[1]], organism="human", ont="BP", pvalueCutoff=0.01)
##' 	#head(summary(yy))
##' 	#plot(yy)
##'
enrichGO <- function(gene,
                     organism="human",
                     ont="MF",
                     pvalueCutoff=0.05,
                     qvalueCutoff=0.05,
                     readable=FALSE) {

    enrich.internal(gene,
                    organism=organism,
                    pvalueCutoff=pvalueCutoff,
                    qvalueCutoff=qvalueCutoff,
                    ont=ont,
                    readable=readable)
}


##' @importFrom DOSE EXTID2TERMID
##' @S3method EXTID2TERMID MF
EXTID2TERMID.MF <- function(gene, organism) {
    EXTID2TERMID.GO(gene=gene, ont="MF", organism=organism)
}

##' @importFrom DOSE EXTID2TERMID
##' @S3method EXTID2TERMID BP
EXTID2TERMID.BP <- function(gene, organism) {
    EXTID2TERMID.GO(gene=gene, ont="BP", organism=organism)
}

##' @importFrom DOSE EXTID2TERMID
##' @S3method EXTID2TERMID CC
EXTID2TERMID.CC <- function(gene, organism) {
    EXTID2TERMID.GO(gene=gene, ont="CC", organism=organism)
}

##' @importMethodsFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importFrom plyr dlply
##' @importFrom plyr .
##' @importClassesFrom methods data.frame
EXTID2TERMID.GO <- function(gene, ont, organism) {
    # get all goterms within the specific ontology
    goterms <- Ontology(GOTERM)
    goterms <- names(goterms[goterms == ont])

    ## get organism specific GO terms
    annoDb <- switch(organism,
                     human = "org.Hs.eg.db",
                     mouse = "org.Mm.eg.db",
                     yeast = "org.Sc.sgd.db",
                     )
    require(annoDb, character.only = TRUE)

    mappedDb <- switch(organism,
                     human = "org.Hs.egGO2ALLEGS",
                     mouse = "org.Mm.egGO2ALLEGS",
                     yeast = "org.Sc.sgdGO2ALLORFS",
                     )
    mappedDb <- eval(parse(text=mappedDb))

    orgTerm <- mappedkeys(mappedDb)

    ## narrow down goterms to specific organism
    Terms <- goterms[goterms %in% orgTerm]

    ## mapping GO to External gene ID
    class(Terms) <- ont
    GO2ExtID <- TERMID2EXTID(Terms, organism)


    gene <- as.character(gene)
    qGO2ExtID = lapply(GO2ExtID, function(i) gene[gene %in% i])
    len <- sapply(qGO2ExtID, length)
    notZero.idx <- len != 0
    qGO2ExtID <- qGO2ExtID[notZero.idx]

    len <- sapply(qGO2ExtID, length)
    qGO2ExtID.df <- data.frame(GO=rep(names(qGO2ExtID), times=len),
                               ExtID=unlist(qGO2ExtID))

    ExtID <- NULL ## to satisfy codetools
    qExtID2GO <- dlply(qGO2ExtID.df, .(ExtID), function(i) as.character(i$GO))

    return(qExtID2GO)
}

##' @importFrom DOSE TERMID2EXTID
##' @S3method TERMID2EXTID MF
TERMID2EXTID.MF <- function(term, organism) {
    TERMID2EXTID.GO(term, organism)
}

##' @importFrom DOSE TERMID2EXTID
##' @S3method TERMID2EXTID BP
TERMID2EXTID.BP <- function(term, organism) {
    TERMID2EXTID.GO(term, organism)
}

##' @importFrom DOSE TERMID2EXTID
##' @S3method TERMID2EXTID CC
TERMID2EXTID.CC <- function(term, organism) {
    TERMID2EXTID.GO(term, organism)
}

##' @importMethodsFrom AnnotationDbi mget
TERMID2EXTID.GO <- function(term, organism) {
    term <- as.character(term)
    annoDb <- switch(organism,
                     human = "org.Hs.eg.db",
                     mouse = "org.Mm.eg.db",
                     yeast = "org.Sc.sgd.db",
                     )
    require(annoDb, character.only = TRUE)

    mappedDb <- switch(organism,
                     human = "org.Hs.egGO2ALLEGS",
                     mouse = "org.Mm.egGO2ALLEGS",
                     yeast = "org.Sc.sgdGO2ALLORFS",
                     )
    mappedDb <- eval(parse(text=mappedDb))
    GO2ExtID <- mget(term, mappedDb, ifnotfound=NA)
    GO2ExtID <- lapply(GO2ExtID, function(i) unique(i))
    return(GO2ExtID)
}



##' @importFrom DOSE ALLEXTID
##' @S3method ALLEXTID MF
ALLEXTID.MF <- function(organism) {
    ALLEXTID.GO(organism)
}

##' @importFrom DOSE ALLEXTID
##' @S3method ALLEXTID BP
ALLEXTID.BP <- function(organism) {
    ALLEXTID.GO(organism)
}

##' @importFrom DOSE ALLEXTID
##' @S3method ALLEXTID CC
ALLEXTID.CC <- function(organism) {
    ALLEXTID.GO(organism)
}

##' @importMethodsFrom AnnotationDbi mappedkeys
ALLEXTID.GO <- function(organism) {
    annoDb <- switch(organism,
                     human = "org.Hs.eg.db",
                     mouse = "org.Mm.eg.db",
                     yeast = "org.Sc.sgd.db",
                     )
    require(annoDb, character.only = TRUE)

    mappedDb <- switch(organism,
                       human = "org.Hs.egGO",
                       mouse = "org.Mm.egGO",
                       yeast = "org.Sc.sgdGO",
                       )
    mappedDb <- eval(parse(text=mappedDb))
    extID <- mappedkeys(mappedDb)
    return(extID)
}

##' @importFrom DOSE TERM2NAME
##' @S3method TERM2NAME MF
TERM2NAME.MF <- function(term) {
    TERM2NAME.GO(term)
}

##' @importFrom DOSE TERM2NAME
##' @S3method TERM2NAME BP
TERM2NAME.BP <- function(term) {
    TERM2NAME.GO(term)
}

##' @importFrom DOSE TERM2NAME
##' @S3method TERM2NAME CC
TERM2NAME.CC <- function(term) {
    TERM2NAME.GO(term)
}

##' @importFrom GO.db GOTERM
##' @importMethodsFrom AnnotationDbi Term
##' @importMethodsFrom AnnotationDbi mget
TERM2NAME.GO <- function(term) {
    term <- as.character(term)
    go <- mget(term, GOTERM, ifnotfound=NA)
    termName <- sapply(go, Term)
    return(term)
}
