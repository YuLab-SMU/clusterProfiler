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
##' @param use.KEGG.db whether use KEGG.db for annotation
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
                       organism      ="human",
                       pvalueCutoff  = 0.05,
                       pAdjustMethod ="BH",
                       universe,
                       minGSSize     = 5,
                       qvalueCutoff  =0.2,
                       readable      =FALSE,
                       use.KEGG.db   = TRUE) {

    enrich.internal(gene,
                    organism      = organism,
                    pvalueCutoff  =pvalueCutoff,
                    pAdjustMethod =pAdjustMethod,
                    ont           = "KEGG",
                    universe      = universe,
                    minGSSize     = minGSSize,
                    qvalueCutoff  = qvalueCutoff,
                    readable      = readable,
                    use.KEGG.db)

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
EXTID2TERMID.KEGG <- function(gene, organism, use.KEGG.db=TRUE) {
    gene <- as.character(gene)
    
    if (use.KEGG.db && organism %in% KEGG_supported_organism() ) {
        qExtID2PathID <- mget(gene, KEGGEXTID2PATHID, ifnotfound=NA)
        ## notNA.idx <- unlist(lapply(qExtID2PathID, function(i) !all(is.na(i))))
        ## qExtID2PathID <- qExtID2PathID[notNA.idx]
    } else {
        buildKEGGmap_supported_organism(organism)
        
        if (file.exists("EXTID2KEGGPATHID.rda")) {
            EXTID2KEGGPATHID <- NULL
            load("EXTID2KEGGPATHID.rda")
            qExtID2PathID <- EXTID2KEGGPATHID[gene]
        } else {
            stop("KEGG mapping file not found in the working directory")
        }
        
    }
    removeEmptyEntry.list(qExtID2PathID)
}

##' @importFrom DOSE TERMID2EXTID
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom KEGG.db KEGGPATHID2EXTID
##' @method TERMID2EXTID KEGG
##' @export
TERMID2EXTID.KEGG <- function(term, organism, use.KEGG.db=TRUE) {
    if(use.KEGG.db && organism %in% KEGG_supported_organism()) {
        pathID2ExtID <- mget(unique(term), KEGGPATHID2EXTID, ifnotfound=NA)
    } else {
        buildKEGGmap_supported_organism(organism)
        if(file.exists("KEGGPATHID2EXTID.rda")) {
            KEGGPATHID2EXTID <- NULL
            load("KEGGPATHID2EXTID.rda")
            pathID2ExtID <- KEGGPATHID2EXTID[unique(term)]
        } else {
            stop("KEGG mapping file not found in the working directory")
        }
    }
    removeEmptyEntry.list(pathID2ExtID)
}

##' @importFrom DOSE ALLEXTID
##' @importFrom KEGG.db KEGGPATHID2EXTID
##' @method ALLEXTID KEGG
##' @export
ALLEXTID.KEGG <- function(organism, use.KEGG.db=TRUE) {
    if (use.KEGG.db && organism %in% KEGG_supported_organism()) {
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
        } else if (organism == "worm" || organism == "celegans") {
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
        orgPath2ExtID <- TERMID2EXTID.KEGG(orgPathID, organism, use.KEGG.db)
        orgPath2ExtID <- lapply(orgPath2ExtID, function(i) unique(i))
    } else {
        buildKEGGmap_supported_organism(organism)
        if (file.exists("KEGGPATHID2EXTID.rda")) {
            KEGGPATHID2EXTID <- NULL
            load("KEGGPATHID2EXTID.rda")
            orgPath2ExtID <- KEGGPATHID2EXTID
        } else {
            stop("KEGG mapping file not found in the working directory")
        }
    }
    orgExtID <- unique(unlist(orgPath2ExtID))
    return(orgExtID)
}

##' @importFrom DOSE TERM2NAME
##' @importFrom KEGG.db KEGGPATHID2NAME
##' @importMethodsFrom AnnotationDbi mget
##' @method TERM2NAME KEGG
##' @export
TERM2NAME.KEGG <- function(term, organism, use.KEGG.db=TRUE) {
    term <- as.character(term)

    pathIDs <- gsub("^\\D+", "",term, perl=T)

    if (use.KEGG.db && organism %in% KEGG_supported_organism()) {
        path2name <- unlist(mget(pathIDs, KEGGPATHID2NAME, ifnotfound = NA))
    } else {
        buildKEGGmap_supported_organism(organism)
        if (file.exists("KEGGPATHID2NAME.rda")) {
            KEGGPATHID2NAME <- NULL
            load("KEGGPATHID2NAME.rda")
            path2name <- KEGGPATHID2NAME[pathIDs]
        } else {
            path2name <- TERM2NAME.KEGG(term, organism, TRUE)
        }
    }
    path2name <- path2name[!is.na(path2name)]
    return(path2name)
}

##' download the latest version of KEGG pathway
##'
##' 
##' @title download.KEGG
##' @param species species
##' @return list
##' @author Guangchuang Yu
##' @importFrom KEGGREST keggLink
##' @importFrom KEGGREST keggList
##' @importFrom magrittr %<>%
##' @export
download.KEGG <- function(species) {
    keggpathid2extid <- keggLink(species,"pathway")
    keggpathid2extid %<>% gsub("[^:]+:", "", .)
    names(keggpathid2extid) %<>% gsub("[^:]+:", "", .)

    keggpath2extid.df <- data.frame(pathID=names(keggpathid2extid), extID=keggpathid2extid)
    
    
    ##keggpathway2gene<-split(as.character(keggpathway2gene),names(keggpathway2gene))
    keggpathid2name<-keggList("pathway")
    names(keggpathid2name) %<>% gsub("path:map", "", .)

    res <- list(keggpath2extid = keggpath2extid.df,
                keggpathid2name = keggpathid2name)
    return(res)
}


##' build KEGG annotation files
##'
##' 
##' @title buildKEGGmap
##' @param keggmap pathway to external ID
##' @param id2name pathway id to pathway name
##' @return NULL
##' @importFrom magrittr %>%
##' @author Guangchuang Yu
##' @export
buildKEGGmap <- function(keggmap, id2name=NULL) {
    if (is.null(id2name)) {
        pathid <- keys(KEGGPATHID2NAME)
        id <- keggmap[,1] %>% as.character %>% gsub("^[a-z]+", "", .)
        keggmap <- keggmap[id %in% pathid, ]
    } else {
        KEGGPATHID2NAME <- id2name
        save(KEGGPATHID2NAME, file="KEGGPATHID2NAME.rda")
    }
    
    KEGGPATHID2EXTID <- split(as.character(keggmap[,2]), as.character(keggmap[,1]))
    EXTID2KEGGPATHID <- split(as.character(keggmap[,1]), as.character(keggmap[,2]))
    save(KEGGPATHID2EXTID, file="KEGGPATHID2EXTID.rda")
    save(EXTID2KEGGPATHID, file="EXTID2KEGGPATHID.rda")
}


buildKEGGmap_supported_organism <- function(organism) {
    if (organism == "anopheles") {
        species <- "aga"
    } else if (organism == "arabidopsis") {
        species <- "ath"
    } else if (organism == "bovine") {
        species <- "bta"
    } else if (organism == "canine") {
        species <- "cfa"
    } else if (organism == "chicken") {
        species <- "gga"
    } else if (organism == "chipm") {
        species <- "ptr"
    } else if (organism == "ecolik12") {
        species <- "eco"
    } else if (organism == "ecsakai") {
        species <- "ecs"
    } else if (organism == "fly") {
        species <- "dme"
    } else if (organism == "human") {
        species <- "hsa"
    } else if (organism == "malaria") {
        species <- "pfa"
    } else if (organism == "mouse") {
        species <- "mmu"
    } else if (organism == "pig") {
        species <- "ssc"
    } else if (organism == "rat") {
        species <- "rno"
    } else if (organism == "rhesus") {
        species <- "mcc"
    } else if (organism == "worm" || organism == "celegans") {
        species <- "cel"
    } else if (organism == "xenopus") {
        species <- "xla"
    } else if (organism == "yeast") {
        species <- "sce"
    } else if (organism == "zebrafish") {
        species <- "dre"
    } else {
        species <- NA
    }

    if (!is.na(species)) {
        if (!file.exists("KEGGPATHID2NAME.rda") &&
            !file.exists("KEGGPATHID2EXTID.rda") &&
            !file.exists("EXTID2KEGGPATHID.rda") ) {
            kegg <- download.KEGG(species)
            buildKEGGmap(kegg[[1]], kegg[[2]])
        }
    }
}


KEGG_supported_organism <- function() {
    res <- c(
        "anopheles",
        "arabidopsis",
        "bovine",
        "canine",
        "chicken",
        "chipm",
        "ecolik12",
        "ecsakai",
        "ecs",
        "fly",
        "dme",
        "human",
        "malaria",
        "mouse",
        "pig",
        "rat",
        "rhesus",
        "worm",
        "celegans",
        "xenopus",
        "yeast",
        "zebrafish")
    return(res)
}
