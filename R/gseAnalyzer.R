##' Gene Set Enrichment Analysis of Gene Ontology
##'
##'
##' @title gseGO
##' @param geneList order ranked geneList
##' @param ont one of "BP", "MF", "CC" or "GO"
##' @param organism organism
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


##' Gene Set Enrichment Analysis of KEGG
##'
##'
##' @title gseKEGG
##' @param geneList order ranked geneList
##' @param organism organism
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
gseKEGG <- function(geneList,
                  organism      = "human",
                  exponent      = 1,
                  nPerm         = 1000,
                  minGSSize     = 10,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  verbose       = TRUE) {
    
    gseAnalyzer(geneList      = geneList,
                setType       = "KEGG",
                organism      = organism,
                exponent      = exponent,
                nPerm         = nPerm,
                minGSSize     = minGSSize,
                pvalueCutoff  = pvalueCutoff,
                pAdjustMethod = pAdjustMethod,
                verbose       = verbose)

}

##' visualize analyzing result of GSEA
##'
##' plotting function for gseaResult
##' @title gseaplot
##' @param x gseaResult object
##' @param ... additional parameters
##' @return figure
##' @export
##' @author ygc
gseaplot <- function(x, ...) {
    plot(x, type="gseaplot", ...)
}

##' @title getGeneSet.KEGG
##' @param setType gene set type
##' @param organism organism
##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi as.list
##' @method getGeneSet KEGG
##' @export
getGeneSet.KEGG <- function(setType="KEGG", organism) {
    gs <- as.list(KEGGPATHID2EXTID)
    return(gs)
}


##' @title getGeneSet.GO
##' @param setType gene set type
##' @param organism organism
##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi as.list
##' @method getGeneSet GO
##' @export
getGeneSet.GO <- function(setType="GO", organism) {
    GO2ALLEG <- GO2EXTID(organism)
    gs <- as.list(GO2ALLEG)
    return(gs)
}

##' @title getGeneSet.BP
##' @param setType gene set type
##' @param organism organism
##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi mget
##' @importFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @method getGeneSet BP
##' @export
getGeneSet.BP <- function(setType="BP", organism) {
    gs <- getGeneSet.GO("GO", organism)
    ont <- lapply(mget(names(gs), GOTERM), Ontology)
    ont <- unlist(ont)
    gs <- gs[ont == setType]
    return(gs)
}

##' @title getGeneSet.MF
##' @param setType gene set type
##' @param organism organism
##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi mget
##' @importFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @method getGeneSet MF
##' @export
getGeneSet.MF <- function(setType="MF", organism) {
    gs <- getGeneSet.GO("GO", organism)
    ont <- lapply(mget(names(gs), GOTERM), Ontology)
    ont <- unlist(ont)
    gs <- gs[ont == setType]
    return(gs)
}

##' @title getGeneSet.CC
##' @param setType gene set type
##' @param organism organism
##' @importFrom DOSE getGeneSet
##' @importFrom AnnotationDbi mget
##' @importFrom AnnotationDbi Ontology
##' @importFrom GO.db GOTERM
##' @method getGeneSet CC
##' @export
getGeneSet.CC <- function(setType="CC", organism) {
    gs <- getGeneSet.GO("GO", organism)
    ont <- lapply(mget(names(gs), GOTERM), Ontology)
    ont <- unlist(ont)
    gs <- gs[ont == setType]
    return(gs)
}
