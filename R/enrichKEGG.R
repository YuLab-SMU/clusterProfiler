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
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param qvalueCutoff qvalue cutoff
##' @param readable whether mapping gene ID to gene Name
##' @param use_internal_data logical, if TRUE, use KEGG.db.
##'                default is FALSE, will download online KEGG data
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
                       organism          = "human",
                       pvalueCutoff      = 0.05,
                       pAdjustMethod     = "BH",
                       universe,
                       minGSSize         = 5,
                       qvalueCutoff      = 0.2,
                       readable          = FALSE,
                       use_internal_data = FALSE) {

    enrich.internal(gene,
                    organism      = organism,
                    pvalueCutoff  =pvalueCutoff,
                    pAdjustMethod =pAdjustMethod,
                    ont           = "KEGG",
                    universe      = universe,
                    minGSSize     = minGSSize,
                    qvalueCutoff  = qvalueCutoff,
                    readable      = readable,
                    use_internal_data)

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
## @importFrom KEGG.db KEGGEXTID2PATHID
##' @method EXTID2TERMID KEGG
##' @export
EXTID2TERMID.KEGG <- function(gene, organism, ...) {
    EXTID2TERMID.KEGG.internal(gene, organism, ...)
}

EXTID2TERMID.KEGG.internal <- function(gene, organism, use_internal_data=TRUE, ...) {
    gene <- as.character(gene)
    organism <- organismMapper(organism)
    
    if (use_internal_data && organism %in% KEGG_db_supported() ) {
        KEGGEXTID2PATHID <- get_KEGG_db("KEGGEXTID2PATHID")
        qExtID2PathID <- mget(gene, KEGGEXTID2PATHID, ifnotfound=NA)
     } else {
        EXTID2KEGGPATHID <- get_KEGG_Anno(organism, "EXTID2KEGGPATHID")
        qExtID2PathID <- EXTID2KEGGPATHID[gene]
     }
    removeEmptyEntry.list(qExtID2PathID)
}

##' @importFrom DOSE TERMID2EXTID
##' @importMethodsFrom AnnotationDbi mget
## @importFrom KEGG.db KEGGPATHID2EXTID
##' @method TERMID2EXTID KEGG
##' @export
TERMID2EXTID.KEGG <- function(term, organism, ...) {
    TERMID2EXTID.KEGG.internal(term, organism, ...)
}

TERMID2EXTID.KEGG <- function(term, organism, use_internal_data=TRUE, ...) {
    organism <- organismMapper(organism)
    if(use_internal_data && organism %in% KEGG_db_supported()) {
        KEGGPATHID2EXTID <- get_KEGG_db("KEGGPATHID2EXTID")
        pathID2ExtID <- mget(unique(term), KEGGPATHID2EXTID, ifnotfound=NA)
    } else {
        KEGGPATHID2EXTID <- get_KEGG_Anno(organism, "KEGGPATHID2EXTID")
        pathID2ExtID <- KEGGPATHID2EXTID[unique(term)]
    }
    removeEmptyEntry.list(pathID2ExtID)
}

##' @importFrom DOSE ALLEXTID
## @importFrom KEGG.db KEGGPATHID2EXTID
##' @method ALLEXTID KEGG
##' @export
ALLEXTID.KEGG <- function(organism, ...) {
    ALLEXTID.KEGG.internal(organism, ...)
}

ALLEXTID.KEGG.internal <- function(organism, use_internal_data=TRUE, ...) {
    organism <- organismMapper(organism)
    
    if (use_internal_data && organism %in% KEGG_db_supported()) {
        KEGGPATHID2EXTID <- get_KEGG_db("KEGGPATHID2EXTID")
        pathID <- mappedkeys(KEGGPATHID2EXTID)
        idx <- grep(paste0("^", organism), pathID)
        orgPathID <- pathID[idx]
        class(orgPathID) <- "KEGG"
        orgPath2ExtID <- TERMID2EXTID.KEGG(orgPathID, organism, use_internal_data)
        orgPath2ExtID <- lapply(orgPath2ExtID, function(i) unique(i))
    } else {
        orgPath2ExtID <- get_KEGG_Anno(organism, "KEGGPATHID2EXTID")
    }
    orgExtID <- unique(unlist(orgPath2ExtID))
    return(orgExtID)
}

##' @importFrom DOSE TERM2NAME
## @importFrom KEGG.db KEGGPATHID2NAME
##' @importMethodsFrom AnnotationDbi mget
##' @method TERM2NAME KEGG
##' @export
TERM2NAME.KEGG <- function(term, organism, ...) {
    TERM2NAME.KEGG.internal(term, organism, ...)
}

TERM2NAME.KEGG.internal <- function(term, organism, use_internal_data=TRUE, ...) {
    term <- as.character(term)
    organism <- organismMapper(organism)
    pathIDs <- gsub("^\\D+", "",term, perl=T)

    if (use_internal_data && organism %in% KEGG_db_supported()) {
        KEGGPATHID2NAME <- get_KEGG_db("KEGGPATHID2NAME")
        path2name <- unlist(mget(pathIDs, KEGGPATHID2NAME, ifnotfound = NA))
    } else {
        KEGGPATHID2NAME <- get_KEGG_Anno(organism, "KEGGPATHID2NAME")
        path2name <- KEGGPATHID2NAME[pathIDs]
    }
    path2name <- path2name[!is.na(path2name)]
    return(path2name)
}



