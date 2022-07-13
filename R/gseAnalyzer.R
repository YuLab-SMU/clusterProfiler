##' Gene Set Enrichment Analysis of Gene Ontology
##'
##'
##' @title gseGO
##' @param geneList order ranked geneList
##' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
##' @param OrgDb an OrgDb or a GSON object.
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
    if (inherits(OrgDb, "GSON")) {
        GO_DATA <- OrgDb
        species <- GO_DATA@species
    } else {
        GO_DATA <- get_GO_data(OrgDb, ont, keyType)
        species <- get_organism(OrgDb)
    }
    

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
    res@organism <- species
    res@setType <- ont
    res@keytype <- keyType

    if (ont == "ALL") {
        if (!inherits(OrgDb, "GSON")){
            res <- add_GO_Ontology(res, GO_DATA)
        } else {
            # do nothing for now
        }    
    }

    return(res)
}



##' Gene Set Enrichment Analysis of KEGG Module
##'
##'
##' @title gseMKEGG
##' @param geneList order ranked geneList
##' @param organism a supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html',
##' or a GSON object.
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
##' @examples
##' \dontrun{
##'   data(geneList, package='DOSE')
##'   mkk <- gseMKEGG(geneList = geneList,
##'                    organism = 'hsa',
##'                    pvalueCutoff = 1)
##'   head(mkk)
##'   kk2 <- gson_KEGG('hsa', KEGG_Type="MKEGG")
##'   mkk2 <- gseMKEGG(geneList = geneList,
##'                   organism = kk2,
##'                   pvalueCutoff = 1)
##'   head(mkk2)
##' }
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

    if (inherits(organism, "character")) {                       
        species <- organismMapper(organism)
        KEGG_DATA <- prepare_KEGG(species, "MKEGG", keyType)
    } else if (inherits(organism, "GSON")) {
        KEGG_DATA <- organism
        species <- KEGG_DATA@species
    } else {
        stop("organism should be a species name or a GSON object")
    }


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
##' @examples
##' \dontrun{
##'   data(geneList, package='DOSE')
##'   gsekegg <- gseKEGG(geneList)
##'   head(gsekegg)
##'   kk <- gson_KEGG('hsa')
##'   gsekegg2 <- gseKEGG(de, organism = kk)
##'   head(gsekegg2)
##' }
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


    if (inherits(organism, "character")) {                       
        species <- organismMapper(organism)
        if (use_internal_data) {
            KEGG_DATA <- get_data_from_KEGG_db(species)
        } else {
            KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
        }
    } else if (inherits(organism, "GSON")) {
        KEGG_DATA <- organism
        species <- KEGG_DATA@species
    } else {
        stop("organism should be a species name or a GSON object")
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



