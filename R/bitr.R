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


bitr_kegg <- function(geneID, fromType, toType, species, drop=TRUE) {
    id_types <- c("ncbi-proteinid", "ncbi-geneid", "uniprot", "kegg")
    fromType <- match.arg(fromType, id_types)
    toType <- match.arg(toType, id_types)

    if (fromType == toType)
        stop("fromType and toType should not be identical...")
    
    if (fromType == "kegg") {
        turl <- paste("http://rest.kegg.jp/conv", toType, species, sep='/')
        tidconv <- kegg_rest(turl)
        idconv <- tidconv
    }

    if (toType == "kegg") {
        furl <- paste("http://rest.kegg.jp/conv", fromType, species, sep='/')
        fidconv <- kegg_rest(furl)
        idconv <- fidconv
    }

    if (fromType != "kegg" && toType != "kegg") {
        idconv <- merge(fidconv, tidconv, by.x='from', by.y='from')
        idconv <- idconv[, -1]
        colnames(idconv) <- c("from", "to")
    }
    
    idconv[,1] %<>% gsub("[^:]+:", "", .)
    idconv[,2] %<>% gsub("[^:]+:", "", .)
    
    res <- idconv[idconv[,1] %in% geneID, ]
    n <- sum(!geneID %in% res[,1])
    if (n > 0) {
        warning(paste0(round(n/length(geneID)*100, 2), "%"), " of input gene IDs are fail to map...")
    }
    
    if (! drop && n > 0) {
        misHit <- data.frame(from = geneID[!geneID %in% res[,1]],
                             to = NA)
        res <- rbind(res, misHit)
    }
    colnames(res) <- c(fromType, toType)
    rownames(res) <- NULL
    return(res)
}
