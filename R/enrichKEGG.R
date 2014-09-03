##' KEGG Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment KEGG
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##'"chicken", "chimp", "ecolik12","ecsakai", "fly", "human",
##'"malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus",
##' "yeast" and "zebrafish".
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance.
##' @export
##' @importFrom DOSE enrich.internal
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @importFrom DOSE setReadable
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
                       pAdjustMethod="BH",
                       universe,
                       minGSSize = 5,
                       qvalueCutoff=0.2,
                       readable=FALSE) {

    enrich.internal(gene,
                    organism = organism,
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = "KEGG",
                    universe = universe,
                    minGSSize = minGSSize,
                    qvalueCutoff = qvalueCutoff,
                    readable = readable)

}

##' viewKEGG function is for visualize KEGG pathways
##' works with enrichResult object to visualize enriched KEGG pathway
##'
##'
##' @param obj enrichResult object
##' @param pathwayID pathway ID or index
##' @param foldChange fold change values
##' @param color.low color of low foldChange genes
##' @param color.high color of high foldChange genes
##' @param kegg.native logical
##' @param out.suffix suffix of output file
## @importFrom pathview pathview
## @importFrom pathview kegg.species.code
##' @references Luo et al. (2013) Pathview: an R/Bioconductor package for 
##'pathway-based data integration and visualization. \emph{Bioinformatics} (Oxford,
##'England), 29:14 1830--1831, 2013. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/content/abstract/29/14/1830.abstract}
##'PMID: 23740750
viewKEGG <- function(obj, pathwayID, foldChange,
                       color.low="green",
                       color.high="red",
                       kegg.native=TRUE,
                       out.suffix="clusterProfiler") {

    if (class(obj) != "enrichResult")
        stop("only enrichResult object supported.")
    if (obj@ontology != "KEGG")
        stop("only KEGG supported.")

    print("viewKEGG is a wrapper function of pathview")
    citation("pathview")
    
    pkg <- "pathview"
    suppressMessages(require(pkg, character.only=TRUE))
    if (is.numeric(pathwayID)) {
        pathwayID <- summary(obj)[pathwayID, 1]
    }
    if (length(pathwayID) == 1 & pathwayID == "all") {
        pathwayID <- summary(obj)[, 1]
    }
    m.fc <- max(abs(foldChange))
    bins <- ceiling(m.fc) * 2
    if (bins < 10)
        bins <- 10
    pathview <- eval(parse(text=pkg))
    res <- lapply(pathwayID, function(pid) {
        pathview(gene.data=foldChange,
                 pathway.id = pid,
                 species = "hsa",
                 limit = list(gene=m.fc, cpd=1),
                 bins = list(gene=bins, cpd=10),
                 low = list(gene=color.low, cpd="blue"),
                 high = list(gene=color.high, cpd="yellow"),
                 kegg.native=kegg.native,
                 out.suffix=out.suffix,
                 new.signature=FALSE)
    })
    return (res)
}

##' @importFrom DOSE EXTID2TERMID
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom KEGG.db KEGGEXTID2PATHID
##' @method EXTID2TERMID KEGG
##' @export
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
##' @method TERMID2EXTID KEGG
##' @export
TERMID2EXTID.KEGG <- function(term, organism) {
    pathID2ExtID <- mget(unique(term), KEGGPATHID2EXTID, ifnotfound=NA)
    return(pathID2ExtID)
}

##' @importFrom DOSE ALLEXTID
##' @importFrom KEGG.db KEGGPATHID2EXTID
##' @method ALLEXTID KEGG
##' @export
ALLEXTID.KEGG <- function(organism) {
    ##pathID2ExtID <- as.list(KEGGPATHID2EXTID)
    ##pathID <- names(pathID2ExtID)
    pathID <- mappedkeys(KEGGPATHID2EXTID)

    ## select species specific pathways
    if (organism == "anopheles") {
        idx <- grep("^aga", pathID)
    } else if (organism == "arabidopsis") {
        idx <- grep("^ath", pathID)
    } else if (organism == "bovine") {
        idx <- grep("^bta", pathID)
    } else if (organism == "canine") {
        idx <- grep("^cfa", pathID)
    } else if (organism == "chicken") {
        idx <- grep("^gga", pathID)
    } else if (organism == "chipm") {
        idx <- grep("^ptr", pathID)
    } else if (organism == "ecolik12") {
        idx <- grep("^eco", pathID)
    } else if (organism == "ecsakai") {
        idx <- grep("^ecs", pathID)
    } else if (organism == "fly") {
        idx <- grep("^dme", pathID)
    } else if (organism == "human") {
        idx <- grep("^hsa", pathID)
    } else if (organism == "malaria") {
        idx <- grep("^pfa", pathID)
    } else if (organism == "mouse") {
        idx <- grep("^mmu", pathID)
    } else if (organism == "pig") {
        idx <- grep("^ssc", pathID)
    } else if (organism == "rat") {
        idx <- grep("^rno", pathID)
    } else if (organism == "rhesus") {
        idx <- grep("^mcc", pathID)
    } else if (organism == "worm") {
        idx <- grep("^cel", pathID)
    } else if (organism == "xenopus") {
        idx <- grep("^xla", pathID)
    } else if (organism == "yeast") {
        idx <- grep("^sce", pathID)
    } else if (organism == "zebrafish") {
        idx <- grep("^dre", pathID)
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
##' @method TERM2NAME KEGG
##' @export
TERM2NAME.KEGG <- function(term, organism) {
    term <- as.character(term)
    pathIDs <- gsub("^\\D+", "",term, perl=T)
    path2name <- unlist(mget(pathIDs, KEGGPATHID2NAME))
    return(path2name)
}

