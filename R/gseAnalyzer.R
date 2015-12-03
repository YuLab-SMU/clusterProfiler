##' Gene Set Enrichment Analysis of Gene Ontology
##'
##'
##' @title gseGO
##' @param geneList order ranked geneList
##' @param ont one of "BP", "MF", "CC" or "GO"
##' @param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##'"chicken", "chimp", "coelicolor", "ecolik12","ecsakai", "fly", "gondii","human",
##'"malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus", "yeast" and
##'"zebrafish".
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @importFrom DOSE gsea
##' @importFrom DOSE gseAnalyzer
##' @importClassesFrom DOSE gseaResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseGO <- function(geneList,
                  ont           = "BP", 
                  organism      = "human",
                  exponent      = 1,
                  nPerm         = 1000,
                  minGSSize     = 10,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  verbose       = TRUE) {
    
    gseAnalyzer(geneList      = geneList,
                setType       = ont,
                organism      = organism,
                exponent      = exponent,
                nPerm         = nPerm,
                minGSSize     = minGSSize,
                pvalueCutoff  = pvalueCutoff,
                pAdjustMethod = pAdjustMethod,
                verbose       = verbose)
    
}

##' Gene Set Enrichment Analysis of KEGG Module
##'
##'
##' @title gseMKEGG
##' @param geneList order ranked geneList
##' @param organism all KEGG Module supported organisms
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseMKEGG <- function(geneList,
                    organism          = "human",
                    exponent          = 1,
                    nPerm             = 1000,
                    minGSSize         = 10,
                    pvalueCutoff      = 0.05,
                    pAdjustMethod     = "BH",
                    verbose           = TRUE) {
    species <- organismMapper(organism)
    keggModule <- download.KEGG_Module(species)
    res <- GSEA(geneList = geneList,
                exponent = exponent,
                nPerm = nPerm,
                minGSSize = minGSSize,
                pvalueCutoff = pvalueCutoff,
                pAdjustMethod = pAdjustMethod,
                TERM2GENE = keggModule$keggmodule2extid,
                TERM2NAME = keggModule$keggmodule2name,
                verbose = verbose)
    res@setType <- "MKEGG"
    return(res)
}


##' Gene Set Enrichment Analysis of KEGG
##'
##'
##' @title gseKEGG
##' @param geneList order ranked geneList
##' @param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##'"chicken", "chimp", "ecolik12","ecsakai", "fly", "human",
##'"malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus",
##' "yeast" and "zebrafish".
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param use_internal_data whether use KEGG.db or not
##' @param verbose print message or not
##' @importFrom DOSE gsea
##' @importFrom DOSE gseAnalyzer
##' @importClassesFrom DOSE gseaResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseKEGG <- function(geneList,
                    organism          = "human",
                    exponent          = 1,
                    nPerm             = 1000,
                    minGSSize         = 10,
                    pvalueCutoff      = 0.05,
                    pAdjustMethod     = "BH",
                    use_internal_data = FALSE,
                    verbose           = TRUE) {
    
    gseAnalyzer(geneList          = geneList,
                setType           = "KEGG",
                organism          = organism,
                exponent          = exponent,
                nPerm             = nPerm,
                minGSSize         = minGSSize,
                pvalueCutoff      = pvalueCutoff,
                pAdjustMethod     = pAdjustMethod,
                use_internal_data = use_internal_data,
                verbose           = verbose)

}

##' visualize analyzing result of GSEA
##'
##' plotting function for gseaResult
##' @title gseaplot
##' @param gseaResult gseaResult object
##' @param geneSetID geneSet ID
##' @param by one of "runningScore" or "position"
##' @return ggplot2 object
##' @export
##' @author ygc
gseaplot <- DOSE::gseaplot


##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi as.list
##' @method getGeneSet KEGG
##' @export
getGeneSet.KEGG <- function(setType="KEGG", organism, ...) {
    getGeneSet.KEGG.internal(setType, organism, ...)
}

getGeneSet.KEGG.internal <- function(setType="KEGG", organism, use_internal_data=TRUE, ...) {
    if (use_internal_data) {
        KEGGPATHID2EXTID <- get_KEGG_db("KEGGPATHID2EXTID")
        gs <- as.list(KEGGPATHID2EXTID)
    } else {
        gs <- get_KEGG_Anno(organism, "KEGGPATHID2EXTID")
    }
    return(gs)
}


##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi as.list
##' @method getGeneSet GO
##' @export
getGeneSet.GO <- function(setType="GO", organism, ...) {
    GO2ALLEG <- GO2EXTID(organism)
    gs <- as.list(GO2ALLEG)
    return(gs)
}

##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi mget
##' @importFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @method getGeneSet BP
##' @export
getGeneSet.BP <- function(setType="BP", organism, ...) {
    getGeneSet.GO_internal(setType, organism, ...)
}


##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi mget
##' @importFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @method getGeneSet MF
##' @export
getGeneSet.MF <- function(setType="MF", organism, ...) {
    getGeneSet.GO_internal(setType, organism, ...)
}


##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi mget
##' @importFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @method getGeneSet CC
##' @export
getGeneSet.CC <- function(setType="CC", organism, ...) {
    getGeneSet.GO_internal(setType, organism, ...)
}


getGeneSet.GO_internal <- function(setType, organism, ...) {
    gs <- getGeneSet.GO("GO", organism)
    term <- mget(names(gs), GOTERM, ifnotfound=NA)
    gs <- gs[!is.na(term)]
    term <- term[!is.na(term)]
    ont <- lapply(term, Ontology)
    ont <- unlist(ont)
    gs <- gs[ont == setType]
    return(gs)
}
