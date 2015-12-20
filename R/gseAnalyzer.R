##' Gene Set Enrichment Analysis of Gene Ontology
##'
##'
##' @title gseGO
##' @param geneList order ranked geneList
##' @param ont one of "BP", "MF", "CC" or "GO"
##' @param OrgDb OrgDb
##' @param keytype keytype of gene
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @importClassesFrom DOSE gseaResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseGO <- function(geneList,
                  ont           = "BP", 
                  OrgDb,
                  keytype       = "ENTREZID",
                  exponent      = 1,
                  nPerm         = 1000,
                  minGSSize     = 10,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  verbose       = TRUE,
                  seed          = FALSE) {

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
    
    GO_DATA <- get_GO_data(OrgDb, ont, keytype)

    res <-  GSEA_internal(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = GO_DATA,
                          seed = seed)
    res@organism <- get_organism(OrgDb)
    res@setType <- ont
    
    if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA)
    } 
    return(res)
}


##' Gene Set Enrichment Analysis of KEGG Module
##'
##'
##' @title gseMKEGG
##' @param geneList order ranked geneList
##' @param species supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseMKEGG <- function(geneList,
                     species          = 'hsa',
                     exponent          = 1,
                     nPerm             = 1000,
                     minGSSize         = 10,
                     pvalueCutoff      = 0.05,
                     pAdjustMethod     = "BH",
                     verbose           = TRUE,
                     seed = FALSE) {

    KEGG_DATA <- download.KEGG(species, "MKEGG")

    res <-  GSEA_internal(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = KEGG_DATA,
                          seed = seed)
    res@organism <- species
    res@setType <- "MKEGG"

    return(res)
}


##' Gene Set Enrichment Analysis of KEGG
##'
##'
##' @title gseKEGG
##' @inheritParams gseMKEGG
##' @importClassesFrom DOSE gseaResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseKEGG <- function(geneList,
                    species           = 'hsa',
                    exponent          = 1,
                    nPerm             = 1000,
                    minGSSize         = 10,
                    pvalueCutoff      = 0.05,
                    pAdjustMethod     = "BH",
                    verbose           = TRUE,
                    seed              = FALSE) {

    
    KEGG_DATA <- download.KEGG(species, "KEGG")

    res <-  GSEA_internal(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = KEGG_DATA,
                          seed = seed)
    res@organism <- species
    res@setType <- "KEGG"

    return(res)
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


