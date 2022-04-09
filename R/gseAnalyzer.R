##' Gene Set Enrichment Analysis of Gene Ontology
##'
##'
##' @title gseGO
##' @param geneList order ranked geneList
##' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
##' @param OrgDb OrgDb
##' @param keyType keytype of gene
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param eps This parameter sets the boundary for calculating the p value.
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @param ... other parameter
##' @importClassesFrom DOSE gseaResult
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseGO <- function(geneList,
                  ont           = "BP",
                  OrgDb,
                  keyType       = "ENTREZID",
                  exponent      = 1,
                  minGSSize     = 10,
                  maxGSSize     = 500,
                  eps           = 1e-10,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  verbose       = TRUE,
                  seed          = FALSE,
                  by            = 'fgsea',
                  ...) {

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

    GO_DATA <- get_GO_data(OrgDb, ont, keyType)

    res <-  GSEA_internal(geneList      = geneList,
                          exponent      = exponent,
                          minGSSize     = minGSSize,
                          maxGSSize     = maxGSSize,
                          eps           = eps,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose       = verbose,
                          USER_DATA     = GO_DATA,
                          seed          = seed,
                          by            = by,
                          ...)
  
    

    if (is.null(res))
        return(res)
        
    if (keyType == 'SYMBOL') {
        res@readable <- TRUE
    }
    res@organism <- get_organism(OrgDb)
    res@setType <- ont
    res@keytype <- keyType

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
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param eps This parameter sets the boundary for calculating the p value.
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @param ... other parameter
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseMKEGG <- function(geneList,
                     organism          = 'hsa',
                     keyType           = 'kegg',
                     exponent          = 1,
                     minGSSize         = 10,
                     maxGSSize         = 500,
                     eps               = 1e-10,
                     pvalueCutoff      = 0.05,
                     pAdjustMethod     = "BH",
                     verbose           = TRUE,
                     seed              = FALSE,
                     by                = 'fgsea',
                     ...) {

    species <- organismMapper(organism)
    KEGG_DATA <- prepare_KEGG(species, "MKEGG", keyType)
    
    res <-  GSEA_internal(geneList       = geneList,
                          exponent       = exponent,
                          minGSSize      = minGSSize,
                          maxGSSize      = maxGSSize,
                          eps            = eps,
                          pvalueCutoff   = pvalueCutoff,
                          pAdjustMethod  = pAdjustMethod,
                          verbose        = verbose,
                          USER_DATA      = KEGG_DATA,
                          seed           = seed,
                          by             = by,
                          ...)
   

    if (is.null(res))
        return(res)


    res@organism <- species
    res@setType <- "MKEGG"
    res@keytype <- "UNKNOWN"

    return(res)
}


##' Gene Set Enrichment Analysis of KEGG
##'
##'
##' @title gseKEGG
##' @inheritParams gseMKEGG
##' @param use_internal_data logical, use KEGG.db or latest online KEGG data
##' @export
##' @return gseaResult object
##' @author Yu Guangchuang
gseKEGG <- function(geneList,
                    organism          = 'hsa',
                    keyType           = 'kegg',
                    exponent          = 1,
                    minGSSize         = 10,
                    maxGSSize         = 500,
                    eps               = 1e-10,
                    pvalueCutoff      = 0.05,
                    pAdjustMethod     = "BH",
                    verbose           = TRUE,
                    use_internal_data = FALSE,
                    seed              = FALSE,
                    by                = 'fgsea',
                    ...) {

    species <- organismMapper(organism)
    if (use_internal_data) {
        KEGG_DATA <- get_data_from_KEGG_db(species)
    } else {
        KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
    }

    res <-  GSEA_internal(geneList         = geneList,
                          exponent         = exponent,
                          minGSSize        = minGSSize,
                          maxGSSize        = maxGSSize,
                          eps              = eps,
                          pvalueCutoff     = pvalueCutoff,
                          pAdjustMethod    = pAdjustMethod,
                          verbose          = verbose,
                          USER_DATA        = KEGG_DATA,
                          seed             = seed,
                          by               = by,
                          ...)
    

    if (is.null(res))
        return(res)

    res@organism <- species
    res@setType <- "KEGG"
    res@keytype <- "UNKNOWN"

    return(res)
}



