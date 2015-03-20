##' list ID types supported by annoDb
##'
##' 
##' @title idType
##' @param annoDb annotation db
##' @return character vector
##' @importFrom AnnotationDbi keytypes
##' @export
##' @author Guangchuang Yu
idType <- function(annoDb = "org.Hs.eg.db") {
    db <- eval(parse(text=annoDb))
    keytypes(db)
}

##' Biological Id TRanslator
##'
##' 
##' @title bitr
##' @param geneID input gene id
##' @param fromType input id type
##' @param toType output id type
##' @param annoDb annotation db
##' @param drop drop NA or not
##' @return data.frame
##' @importFrom magrittr %>%
##' @importFrom magrittr %<>%
##' @importFrom AnnotationDbi select
##' @export
##' @author Guangchuang Yu
bitr <- function(geneID, fromType, toType, annoDb, drop=TRUE) {
    geneID %<>% as.character %>% unique
    db <- eval(parse(text=annoDb))
    res <- suppressWarnings(select(db,
                                   keys = geneID,
                                   keytype = fromType, 
                                   columns=c(fromType, toType)))
    
    ii <- which(is.na(res[,2]))
    n <- res[ii, 1] %>% unique %>% length
    if (n) {
        warning(paste0(n/length(geneID)*100, "%"), "of input gene IDs are fail to map...")
    }
    
    if (drop) {
        res <- res[-ii, ]
    }
    return(res)
}


