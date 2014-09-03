##' GO Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment GO
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##'"chicken", "chimp", "coelicolor", "ecolik12","ecsakai", "fly", "gondii","human",
##'"malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus", "yeast" and
##'"zebrafish".
##' @param ont One of "MF", "BP", and "CC" subontologies.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance.
##' @importFrom DOSE enrich.internal
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @importFrom DOSE setReadable
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
                     pAdjustMethod="BH",
                     universe,
                     qvalueCutoff = 0.2,
                     minGSSize = 5,
                     readable=FALSE) {

    enrich.internal(gene,
                    organism=organism,
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    universe = universe,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize = minGSSize,
                    readable = readable)
}


##' enrichment map
##'
##' enrichMap
##' @title enrichMap
##' @param x enrichResult or gseaResult
##' @param ... additional parameter
##' @return figure
##' @export
##' @author ygc
enrichMap <- function(x, ...) {
    ## DOSE::enrichMap(...)
    plot(x, type="enrichMap", ...)
}

##' category-gene-net plot
##'
##' category gene association
##' @title cnetplot
##' @param x enrichResult object
##' @param ... additional parameter
##' @return figure
##' @export
##' @author ygc
cnetplot <- function(x, ...) {
    plot(x, type="cnet", ...)
}

##' @importFrom DOSE EXTID2TERMID
##' @method EXTID2TERMID MF
##' @export
EXTID2TERMID.MF <- function(gene, organism) {
    EXTID2TERMID.GO(gene=gene, ont="MF", organism=organism)
}

##' @importFrom DOSE EXTID2TERMID
##' @method EXTID2TERMID BP
##' @export
EXTID2TERMID.BP <- function(gene, organism) {
    EXTID2TERMID.GO(gene=gene, ont="BP", organism=organism)
}

##' @importFrom DOSE EXTID2TERMID
##' @method EXTID2TERMID CC
##' @export
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

    gene <- as.character(gene)

    ## get all goterms within the specific ontology
    goterms <- Ontology(GOTERM)
    goterms <- names(goterms[goterms == ont])

    supported_Org <- getSupported_Org()
    if (organism %in% supported_Org) {
        mappedDb <- getGO2ALLEG_MappedDb(organism)

        orgTerm <- mappedkeys(mappedDb)

        ## narrow down goterms to specific organism
        Terms <- goterms[goterms %in% orgTerm]

        ## mapping GO to External gene ID
        class(Terms) <- ont
        GO2ExtID <- TERMID2EXTID(Terms, organism)


        qGO2ExtID = lapply(GO2ExtID, function(i) gene[gene %in% i])
        len <- sapply(qGO2ExtID, length)
        notZero.idx <- len != 0
        qGO2ExtID <- qGO2ExtID[notZero.idx]

        len <- sapply(qGO2ExtID, length)
        qGO2ExtID.df <- data.frame(GO=rep(names(qGO2ExtID), times=len),
                                   ExtID=unlist(qGO2ExtID))

        ExtID <- NULL ## to satisfy codetools
        qExtID2GO <- dlply(qGO2ExtID.df, .(ExtID), function(i) as.character(i$GO))
    } else {
        oldwd <- getwd()
        if(organism == "D39") {
            dir <- system.file("extdata/D39/", package="clusterProfiler")
            setwd(dir)
        }
        if(organism == "M5005") {
            dir <- system.file("extdata/M5005/", package="clusterProfiler")
            setwd(dir)
        }
        if (file.exists("EG2ALLGO.rda")) {
            EG2ALLGO <- NULL # to satisfy codetools
            load("EG2ALLGO.rda")
            qExtID2GO <- EG2ALLGO[gene]
            qExtID2GO <- lapply(qExtID2GO, function(i) i[i %in% goterms])
        } else {
            setwd(oldwd)
            stop("GO mapping file not found in the working directory")
        }
        setwd(oldwd)
    }
    return(qExtID2GO)
}

##' @importFrom DOSE TERMID2EXTID
##' @method TERMID2EXTID MF
##' @export
TERMID2EXTID.MF <- function(term, organism) {
    TERMID2EXTID.GO(term, organism)
}

##' @importFrom DOSE TERMID2EXTID
##' @method TERMID2EXTID BP
##' @export
TERMID2EXTID.BP <- function(term, organism) {
    TERMID2EXTID.GO(term, organism)
}

##' @importFrom DOSE TERMID2EXTID
##' @method TERMID2EXTID CC
##' @export
TERMID2EXTID.CC <- function(term, organism) {
    TERMID2EXTID.GO(term, organism)
}

##' @importMethodsFrom AnnotationDbi mget
##' @importFrom GOSemSim getSupported_Org
TERMID2EXTID.GO <- function(term, organism) {
    term <- as.character(term)

    GO2ALLEG <- GO2EXTID(organism)
    if (is(GO2ALLEG, "Go3AnnDbBimap")) {
        GO2ExtID <- mget(term, GO2ALLEG, ifnotfound=NA)
        GO2ExtID <- lapply(GO2ExtID, function(i) unique(i))
    } else {
        GO2ExtID <- GO2ALLEG[term]
    }
    return(GO2ExtID)
}

GO2EXTID <- function(organism) {
    supported_Org <- getSupported_Org()
    if (organism %in% supported_Org) {
        GO2ALLEG <- getGO2ALLEG_MappedDb(organism)
    } else {
        oldwd <- getwd()
        if(organism == "D39") {
            dir <- system.file("extdata/D39/", package="clusterProfiler")
            setwd(dir)
        }
        if(organism == "M5005") {
            dir <- system.file("extdata/M5005/", package="clusterProfiler")
            setwd(dir)
        }
        if (file.exists("GO2ALLEG.rda")) {
            GO2ALLEG <- NULL # to satisfy codetools
            load("GO2ALLEG.rda")
        } else {
            setwd(oldwd)
            stop("GO Mapping file not found in the working directory")
        }
        setwd(oldwd)
    }
    return(GO2ALLEG)
}


##' @importFrom DOSE ALLEXTID
##' @method ALLEXTID MF
##' @export
ALLEXTID.MF <- function(organism) {
    ALLEXTID.GO(organism)
}

##' @importFrom DOSE ALLEXTID
##' @method ALLEXTID BP
##' @export
ALLEXTID.BP <- function(organism) {
    ALLEXTID.GO(organism)
}

##' @importFrom DOSE ALLEXTID
## @S3method ALLEXTID CC
##' @method ALLEXTID CC
##' @export
ALLEXTID.CC <- function(organism) {
    ALLEXTID.GO(organism)
}


##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importFrom GOSemSim getSupported_Org
ALLEXTID.GO <- function(organism) {
    supported_Org <- getSupported_Org()
    if (organism %in% supported_Org) {
        mappedDb <- getEG2GO_MappedDb(organism)
        extID <- mappedkeys(mappedDb)
    } else {
        oldwd <- getwd()
        if(organism == "D39") {
            dir <- system.file("extdata/D39/", package="clusterProfiler")
            setwd(dir)
        }
        if(organism == "M5005") {
            dir <- system.file("extdata/M5005/", package="clusterProfiler")
            setwd(dir)
        }
        if (file.exists("EG2ALLGO.rda")) {
            EG2ALLGO <- NULL ## to satisfy codetools
            load("EG2ALLGO.rda")
            extID <- names(EG2ALLGO)
        } else {
            setwd(oldwd)
            stop("GO mapping file not found in the working directory")
        }
        setwd(oldwd)
    }
    return(extID)
}

##' @importFrom DOSE TERM2NAME
##' @method TERM2NAME MF
##' @export
TERM2NAME.MF <- function(term, organism) {
    TERM2NAME.GO(term, organism)
}

##' @importFrom DOSE TERM2NAME
##' @method TERM2NAME BP
##' @export
TERM2NAME.BP <- function(term, organism) {
    TERM2NAME.GO(term, organism)
}

##' @importFrom DOSE TERM2NAME
##' @method TERM2NAME CC
##' @export
TERM2NAME.CC <- function(term, organism) {
    TERM2NAME.GO(term, organism)
}

##' @importFrom GO.db GOTERM
##' @importMethodsFrom AnnotationDbi Term
##' @importMethodsFrom AnnotationDbi mget
TERM2NAME.GO <- function(term, organism) {
    term <- as.character(term)
    go <- mget(term, GOTERM, ifnotfound=NA)
    termName <- sapply(go, Term)
    return(termName)
}
