
mydownload <- function(url, method = NULL, quiet = TRUE, ...) {
    if (is.null(method))
        method <- getOption("clusterProfiler.download.method")

    if (!is.null(method) && method != "auto") { 
        dl <- tryCatch(utils::download.file(url, quiet = TRUE, method = method, ...),
                       error = function(e) NULL)       
    } else {
        dl <- tryCatch(downloader::download(url, quiet = TRUE, ...),
                       error = function(e) NULL)        
    }

    return(dl)
}



GI2EG <- function(GI, organism="D39") {
    gi <- as.character(GI)
    ## remove blank
    gi <- sub("^\\s+", "", gi, perl=T)
    gi <- sub("\\s+$", "", gi, perl=T)
    ## remove GI: or gi:
    gi <- sapply(gi, function(i) unlist(strsplit(i, split="\\|"))[2])
    ## load corresponding gene table and protein table
    geneTable <- proteinTable <- NULL
    if (organism == "M5005") {
        f <- system.file("extdata", "M5005/geneTable.rda", package="clusterProfiler")
        load(f)
        gi.eg <- geneTable[geneTable$GI %in% gi, c("GI", "GeneID")]
    } else if (organism=="D39") {
        gt <- system.file("extdata", "D39/geneTable.rda", package="clusterProfiler")
        load(gt)
        pt <- system.file("extdata", "D39/proteinTable.rda", package="clusterProfiler")
        load(pt)
        idx <- match(gi, proteinTable$PID)
        locus <- proteinTable[idx, "Synonym"]
        locus <- as.character(locus)
        idx <- match(locus, geneTable$Locus)
        gene <- geneTable[idx, "GeneID"]
        gene <- as.character(gene)
        gi.eg <- data.frame(GI=gi, GeneID=gene)
    } else {
        stop("not supported yet...")
    }

    return(gi.eg)
}


removeEmptyEntry.list <- function(x) {
    notNA.idx <- unlist(lapply(x, function(i) !is.null(i) && !all(is.na(i))))
    x[notNA.idx]
}


GSEA_internal <- DOSE:::GSEA_internal
enricher_internal <- DOSE:::enricher_internal
get_fun_from_pkg <- yulab.utils:::get_fun_from_pkg

globalVariables(".")
