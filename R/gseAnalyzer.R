#' Gene Set Enrichment Analysis of Gene Ontology
#'
#'
#' @title gseGO
#' @param geneList order ranked geneList
#' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param OrgDb OrgDb
#' @param keyType keytype of gene
#' @param exponent weight of each step
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param nPerm The number of permutations.
#' @param method method of calculating the pvalue, one of "multilevel", "permute" and "sample"
#' @param adaptive logical, whether to use adaptive method for calculating pvalue
#' @param minPerm minimal number of permutations for adaptive method
#' @param maxPerm maximal number of permutations for adaptive method
#' @param pvalThreshold pvalue threshold for adaptive method
#' @param eps boundary for calculating the p value in multilevel mode
#' @param pvalueCutoff pvalue Cutoff
#' @param pAdjustMethod pvalue adjustment method
#' @param verbose print message or not
#' @param seed random seed for reproducibility, set to a number (or TRUE to use a
#'   fixed default seed) to make the result reproducible, or FALSE (default) to draw
#'   a random seed on each run, so results may vary between runs. The underlying
#'   permutation engine uses its own RNG seeded with this value; see
#'   \code{enrichit::gsea()} for details.
#' @param ... other parameters passed to \code{enrichit::gsea_gson()}
#' @importClassesFrom enrichit gseaResult
#' @export
#' @return gseaResult object
#' @author Yu Guangchuang
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
                  nPerm         = 1000,
                  method        = "multilevel",
                  adaptive      = FALSE,
                  minPerm       = 101,
                  maxPerm       = 1e5,
                  pvalThreshold = 0.1,
                  seed          = FALSE,
                  ...) {

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

    GO_DATA <- get_GO_data(OrgDb, ont, keyType)

    res <-  enrichit::gsea_gson(geneList      = geneList,
                          exponent      = exponent,
                          minGSSize     = minGSSize,
                          maxGSSize     = maxGSSize,
                          eps           = eps,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose       = verbose,
                          gson          = GO_DATA,
                          nPerm         = nPerm,
                          method        = method,
                          adaptive      = adaptive,
                          minPerm       = minPerm,
                          maxPerm       = maxPerm,
                          pvalThreshold = pvalThreshold,
                          seed          = seed,
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



#' Gene Set Enrichment Analysis of KEGG Module
#'
#'
#' @title gseMKEGG
#' @param geneList order ranked geneList
#' @param organism supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
#' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
#' @param exponent weight of each step
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param nPerm The number of permutations.
#' @param method method of calculating the pvalue, one of "multilevel", "permute" and "sample"
#' @param adaptive logical, whether to use adaptive method for calculating pvalue
#' @param minPerm minimal number of permutations for adaptive method
#' @param maxPerm maximal number of permutations for adaptive method
#' @param pvalThreshold pvalue threshold for adaptive method
#' @param eps boundary for calculating the p value in multilevel mode
#' @param pvalueCutoff pvalue Cutoff
#' @param pAdjustMethod pvalue adjustment method
#' @param verbose print message or not
#' @param seed random seed for reproducibility, set to a number (or TRUE to use a
#'   fixed default seed) to make the result reproducible, or FALSE (default) to draw
#'   a random seed on each run, so results may vary between runs. The underlying
#'   permutation engine uses its own RNG seeded with this value; see
#'   \code{enrichit::gsea()} for details.
#' @param ... other parameters passed to \code{enrichit::gsea_gson()}
#' @export
#' @return gseaResult object
#' @author Yu Guangchuang
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
                     nPerm             = 1000,
                     method            = "multilevel",
                     adaptive          = FALSE,
                     minPerm           = 101,
                     maxPerm           = 1e5,
                     pvalThreshold     = 0.1,
                     seed              = FALSE,
                     ...) {

    species <- organismMapper(organism)
    KEGG_DATA <- prepare_KEGG(species, "MKEGG", keyType)
    
    res <-  enrichit::gsea_gson(geneList       = geneList,
                          exponent       = exponent,
                          minGSSize      = minGSSize,
                          maxGSSize      = maxGSSize,
                          eps            = eps,
                          pvalueCutoff   = pvalueCutoff,
                          pAdjustMethod  = pAdjustMethod,
                          verbose        = verbose,
                          gson           = KEGG_DATA,
                          nPerm          = nPerm,
                          method         = method,
                          adaptive       = adaptive,
                          minPerm        = minPerm,
                          maxPerm        = maxPerm,
                          pvalThreshold  = pvalThreshold,
                          seed           = seed,
                          ...)
   

    if (is.null(res))
        return(res)


    res@organism <- species
    res@setType <- "MKEGG"
    res@keytype <- "UNKNOWN"

    res <- append_kegg_category(res)
    return(res)
}


#' Gene Set Enrichment Analysis of KEGG
#'
#'
#' @title gseKEGG
#' @inheritParams gseMKEGG
#' @param use_internal_data logical, use KEGG.db or latest online KEGG data
#' @export
#' @return gseaResult object
#' @author Yu Guangchuang
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
                    nPerm             = 1000,
                    method            = "multilevel",
                    adaptive          = FALSE,
                    minPerm           = 101,
                    maxPerm           = 1e5,
                    pvalThreshold     = 0.1,
                    seed              = FALSE,
                    ...) {

    if (inherits(organism, "character")) {           
        if (organism == "cpd") {
            organism = gson_cpd()
        }
    }

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
        keyType <- KEGG_DATA@keytype
    } else {
        stop("organism should be a species name or a GSON object")
    }


    res <-  enrichit::gsea_gson(geneList         = geneList,
                          exponent         = exponent,
                          minGSSize        = minGSSize,
                          maxGSSize        = maxGSSize,
                          eps              = eps,
                          pvalueCutoff     = pvalueCutoff,
                          pAdjustMethod    = pAdjustMethod,
                          verbose          = verbose,
                          gson             = KEGG_DATA,
                          nPerm            = nPerm,
                          method           = method,
                          adaptive         = adaptive,
                          minPerm          = minPerm,
                          maxPerm          = maxPerm,
                          pvalThreshold    = pvalThreshold,
                          seed             = seed,
                          ...)
    

    if (is.null(res))
        return(res)

    res@organism <- species
    res@setType <- "KEGG"
    #res@keytype <- "UNKNOWN"
    res@keytype <- keyType
    
    return(res)
}



