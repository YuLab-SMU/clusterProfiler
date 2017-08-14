##' list ID types supported by annoDb
##'
##'
##' @title idType
##' @param OrgDb annotation db
##' @return character vector
##' @importFrom GOSemSim load_OrgDb
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

##' convert biological ID using KEGG API
##'
##'
##' @title bitr_kegg
##' @param geneID input gene id
##' @param fromType input id type
##' @param toType output id type
##' @param organism supported organism, can be search using search_kegg_organism function
##' @param drop drop NA or not
##' @return data.frame
##' @export
##' @author Guangchuang Yu
bitr_kegg <- function(geneID, fromType, toType, organism, drop=TRUE) {
    id_types <- c("Path", "Module", "ncbi-proteinid", "ncbi-geneid", "uniprot", "kegg")
    fromType <- match.arg(fromType, id_types)
    toType <- match.arg(toType, id_types)

    if (fromType == toType)
        stop("fromType and toType should not be identical...")

    if (fromType == "Path" || fromType == "Module") {
        idconv <- KEGG_path2extid(geneID, organism, fromType, toType)
    } else if (toType == "Path" || toType == "Module") {
        idconv <- KEGG_extid2path(geneID, organism, toType, fromType)
    } else {
        idconv <- KEGG_convert(fromType, toType, organism)
    }

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

KEGG_convert <- function(fromType, toType, species) {
    if (fromType == "kegg" || toType != "kegg") {
        turl <- paste("http://rest.kegg.jp/conv", toType, species, sep='/')
        tidconv <- kegg_rest(turl)
        if (is.null(tidconv))
            stop(toType, " is not supported for ", species, " ...")
        idconv <- tidconv
    }

    if (toType == "kegg" || fromType != "kegg") {
        furl <- paste("http://rest.kegg.jp/conv", fromType, species, sep='/')
        fidconv <- kegg_rest(furl)
        if (is.null(fidconv))
            stop(fromType, " is not supported for ", species, " ...")
        idconv <- fidconv
    }

    if (fromType != "kegg" && toType != "kegg") {
        idconv <- merge(fidconv, tidconv, by.x='from', by.y='from')
        idconv <- idconv[, -1]
    } else if (fromType != "kegg") {
        idconv <- idconv[, c(2,1)]
    }
    colnames(idconv) <- c("from", "to")

    idconv[,1] %<>% gsub("[^:]+:", "", .)
    idconv[,2] %<>% gsub("[^:]+:", "", .)
    return(idconv)
}


##' query all genes in a KEGG pathway or module
##'
##'
##' @title KEGG_path2extid
##' @param keggID KEGG ID, path or module ID
##' @param species species
##' @param keggType one of 'Path' or 'Module'
##' @param keyType KEGG gene type, one of "ncbi-proteinid", "ncbi-geneid", "uniprot", or "kegg"
##' @return extid vector
##' @author guangchuang yu
KEGG_path2extid <- function(keggID, species=sub("\\d+$", "", keggID),
                          keggType = "Path", keyType = "kegg") {
    path2extid <- KEGGPATHID2EXTID(species, keggType, keyType)
    path2extid[path2extid$from %in% keggID, ]
}

KEGG_extid2path <- function(geneID, species, keggType = "Path", keyType = "kegg") {
    path2extid <- KEGGPATHID2EXTID(species, keggType, keyType)
    res <- path2extid[path2extid$to %in% geneID, ]
    res <- res[, c(2,1)]
    colnames(res) <- colnames(path2extid)
    return(res)
}


KEGGPATHID2EXTID <- function(species, keggType = "Path", keyType = "kegg") {
    keggType <- match.arg(keggType, c("Path", "Module"))
    if (keggType == "Path") {
        keggType <- "KEGG"
    } else {
        keggType <- "MKEGG"
    }
    kegg <- download_KEGG(species, keggType, keyType)
    return(kegg$KEGGPATHID2EXTID)
}
