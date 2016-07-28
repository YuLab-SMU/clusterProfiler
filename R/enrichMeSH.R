##' MeSH term enrichment analysis
##'
##' 
##' @title enrichMeSH
##' @param gene a vector of entrez gene id
##' @param MeSHDb MeSHDb
##' @param database one of 'gendoo', 'gene2pubmed' or 'RBBH'
##' @param category one of "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N", "V", "Z"
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @return An \code{enrichResult} instance.
##' @export
##' @author Guangchuang Yu
enrichMeSH <- function(gene,
                       MeSHDb,
                       database = 'gendoo',
                       category = 'C',
                       pvalueCutoff=0.05,
                       pAdjustMethod="BH",
                       universe,
                       qvalueCutoff = 0.2,
                       minGSSize = 10,
                       maxGSSize = 500) {

    MeSH_DATA <- get_MeSH_data(MeSHDb, database, category)
    
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = MeSH_DATA
                             )
    
    meshdb <- get_fun_from_pkg("MeSH.db", "MeSH.db")
    id <- res@result$ID
    mesh2name <- select(meshdb, keys=id, columns=c('MESHID', 'MESHTERM'), keytype='MESHID')
    res@result$Description <- mesh2name[match(id, mesh2name[,1]), 2]
    res@organism <- get_organism(MeSHDb)
    res@ontology <- "MeSH"
    
    return(res)
}


get_MeSH_data <- function(MeSHDb, database, category) {
    category %<>% toupper
    categories <- c("A", "B", "C", "D",
                    "E", "F", "G", "H",
                    "I", "J", "K", "L",
                    "M", "N", "V", "Z")
    
    if (!all(category %in% categories)) {
        stop("please check your 'category' parameter...")
    }
    
    db <- get_fun_from_pkg(MeSHDb, MeSHDb)
    mesh <- select(db, keys=database, columns = c("GENEID", "MESHID","MESHCATEGORY"), keytype = "SOURCEDB")

    mesh <- mesh[ mesh[,3] %in% category, ] 
    mesh2gene <- mesh[, c(2,1)]

    ## meshdb <- get_fun_from_pkg("MeSH.db", "MeSH.db")
    ## mesh2name <- select(meshdb, keys=unique(mesh2gene[,1]), columns=c('MESHID', 'MESHTERM'), keytype='MESHID')

    build_Anno(mesh2gene)
}


get_fun_from_pkg <- function(pkg, fun) {
    ## requireNamespace(pkg)
    ## eval(parse(text=paste0(pkg, "::", fun)))
    require(pkg, character.only = TRUE)
    eval(parse(text = fun))
}


                       
