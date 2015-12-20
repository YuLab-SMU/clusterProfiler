##' KEGG Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment KEGG
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param species supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param qvalueCutoff qvalue cutoff
##' @return A \code{enrichResult} instance.
##' @export
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
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
                       species           = "hsa",
                       pvalueCutoff      = 0.05,
                       pAdjustMethod     = "BH",
                       universe,
                       minGSSize         = 5,
                       qvalueCutoff      = 0.2) {

    KEGG_DATA <- download.KEGG(species, "KEGG")
    res <- enricher_internal(gene,
                             pvalueCutoff  =pvalueCutoff,
                             pAdjustMethod =pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             USER_DATA = KEGG_DATA)
    res@ontology <- "KEGG"
    res@organism <- species
    return(res)
}

get_KEGG_Env <- function() {
    if (! exists("KEGG_clusterProfiler_Env", envir = .GlobalEnv)) {
        assign("KEGG_clusterProfiler_Env", new.env(), .GlobalEnv)
    }
    get("KEGG_clusterProfiler_Env", envir = .GlobalEnv)
}

##' download the latest version of KEGG pathway
##'
##' 
##' @title download.KEGG
##' @param species species
##' @param KEGG_Type one of 'KEGG' or 'MKEGG'
##' @return list
##' @author Guangchuang Yu
##' @importFrom KEGGREST keggLink
##' @importFrom KEGGREST keggList
##' @importFrom magrittr %<>%
download.KEGG <- function(species, KEGG_Type="KEGG") {
    KEGG_Env <- get_KEGG_Env()
    
    use_cached <- FALSE
    
    if (exists("organism", envir = KEGG_Env, inherits = FALSE) &&
        exists("_type_", envir = KEGG_Env, inherits = FALSE) ) {
        
        org <- get("organism", envir=KEGG_Env)
        type <- get("_type_", envir=KEGG_Env)
        
        if (org == species && type == KEGG_Type &&
            exists("KEGGPATHID2NAME", envir=KEGG_Env, inherits = FALSE) &&
            exists("KEGGPATHID2EXTID", envir=KEGG_Env, inherits = FALSE)) {
            
            use_cached <- TRUE
        } 
    }
    
    if (use_cached) {
        KEGGPATHID2EXTID <- get("KEGGPATHID2EXTID", envir=KEGG_Env)
        KEGGPATHID2NAME <- get("KEGGPATHID2NAME", envir=KEGG_Env)
    } else {
        if (KEGG_Type == "KEGG") {
            kres <- download.KEGG.Path(species)
        } else {
            kres <- download.KEGG.Module(species)
        }
        
        KEGGPATHID2EXTID <- kres$KEGGPATHID2EXTID
        KEGGPATHID2NAME <- kres$KEGGPATHID2NAME

        assign("organism", species, envir=KEGG_Env)
        assign("_type_", KEGG_Type, envir=KEGG_Env)
        assign("KEGGPATHID2NAME", KEGGPATHID2NAME, envir=KEGG_Env)
        assign("KEGGPATHID2EXTID", KEGGPATHID2EXTID, envir=KEGG_Env)
    }
        
    build_Anno(KEGGPATHID2EXTID,
               KEGGPATHID2NAME)
}

download.KEGG.Path <- function(species) {
    keggpathid2extid <- tryCatch(keggLink(species,"pathway"), error=function(e) NULL)
    
    if (is.null(keggpathid2extid)) {
        stop("'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...")
    }
    
    keggpathid2extid %<>% gsub("[^:]+:", "", .)
    names(keggpathid2extid) %<>% gsub("[^:]+:", "", .)
    
    keggpath2extid.df <- data.frame(pathID=names(keggpathid2extid), extID=keggpathid2extid)
    
    keggpathid2name <- keggList("pathway")
    names(keggpathid2name) %<>% gsub("path:map", species, .)
    keggpathid2name.df <- data.frame(keggID=names(keggpathid2name),
                                     keggName=keggpathid2name)
    return(list(KEGGPATHID2EXTID=keggpath2extid.df,
                KEGGPATHID2NAME=keggpathid2name.df))
}

download.KEGG.Module <- function(species) {
    keggmodule2extid <- tryCatch(keggLink(species, "module"), error=function(e) NULL)
    if (is.null(keggmodule2extid)) {
        stop("'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...")
    }
    keggmodule2extid %<>% gsub("[^:]+:", "", .)
    names(keggmodule2extid) %<>% gsub("[^:]+:", "", .)
    names(keggmodule2extid) %<>% gsub(species, "", .)
    names(keggmodule2extid) %<>% gsub("^_", "", .)
    keggmodule2extid.df <- data.frame(moduleID=names(keggmodule2extid),
                                      extID = keggmodule2extid)

    keggmodule2name <- keggList("module")
    names(keggmodule2name) %<>% gsub("md:", "", .)
    keggmodule2name.df <- data.frame(moduleID=names(keggmodule2name),
                                     moduleName=keggmodule2name)
    return(list(KEGGPATHID2EXTID=keggmodule2extid.df,
                KEGGPATHID2NAME =keggmodule2name.df))
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






