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
    db <- get_db_obj(annoDb)
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
    idTypes <- idType(annoDb)
    msg <-  paste0("should be one of ", paste(idTypes, collapse=", "), ".")
    if (! fromType %in% idTypes) {
        stop("'fromType' ", msg)
    }
    if (! toType %in% idTypes) {
        stop("'toType' ", msg)
    }
    
    geneID %<>% as.character %>% unique
    db <- get_db_obj(annoDb)
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


get_db_obj <- function(annoDb) {
    require(annoDb, character.only = TRUE)
    eval(parse(text=annoDb))
}

