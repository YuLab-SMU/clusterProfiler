##' list ID types supported by annoDb
##'
##' 
##' @title idType
##' @param OrgDb annotation db
##' @return character vector
##' @importFrom AnnotationDbi keytypes
##' @export
##' @author Guangchuang Yu
idType <- function(OrgDb = "org.Hs.eg.db") {
    db <- load_OrgDb(OrgDb)
    keytypes(db)
}

##' Biological Id TRanslator
##'
##' 
##' @title bitr
##' @param geneID input gene id
##' @param fromType input id type
##' @param toType output id type
##' @param OrgDb annotation db
##' @param drop drop NA or not
##' @return data.frame
##' @importFrom magrittr %>%
##' @importFrom magrittr %<>%
##' @importFrom AnnotationDbi select
##' @export
##' @author Guangchuang Yu
bitr <- function(geneID, fromType, toType, OrgDb, drop=TRUE) {
    idTypes <- idType(OrgDb)
    msg <-  paste0("should be one of ", paste(idTypes, collapse=", "), ".")
    if (! fromType %in% idTypes) {
        stop("'fromType' ", msg)
    }
    if (! all(toType %in% idTypes)) {
        stop("'toType' ", msg)
    }
    
    geneID %<>% as.character %>% unique
    db <- load_OrgDb(OrgDb)
    res <- suppressWarnings(select(db,
                                   keys = geneID,
                                   keytype = fromType, 
                                   columns=c(fromType, toType)))
    
    ii <- which(is.na(res[,2]))
    if (length(ii)) {
        n <- res[ii, 1] %>% unique %>% length
        if (n) {
            warning(paste0(round(n/length(geneID)*100, 2), "%"), " of input gene IDs are fail to map...")
        }
        if (drop) {
            res <- res[-ii, ]
        }
    }
    return(res)
}


