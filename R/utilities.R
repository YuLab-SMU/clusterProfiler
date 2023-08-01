
read_tsv_with_cache <- function(file_url, header = FALSE) {
    read_with_cache(file_url, 
        reader = utils::read.delim, 
        params=list(sep="\t", header = header)
    )
}

# Define a function to read the file with caching
#' @importFrom fs path_join
#' @importFrom digest digest
#' @importFrom memoise memoise
read_with_cache <- memoise::memoise(function(file_url, reader = readLines, params = list()) {
    cache_dir <- getOption("clusterProfiler_cache_dir")
    # Generate a unique cache filename based on the file URL
    cache_filename <- fs::path_join(c(cache_dir, paste0(digest::digest(file_url), ".rds")))
    
    # Check if the cached file exists
    if (file.exists(cache_filename)) {
        # If cached file exists, load and return the cached data
        cached_data <- readRDS(cache_filename)
        return(cached_data)
    } else {
        # If cached file does not exist, download and cache the data
        #data <- read.delim(url(file_url), sep = "\t", header = header)
        data <- do.call(reader, args = c(file_url, params))
        saveRDS(data, cache_filename)
        return(data)
    }
})

# the `read_with_cache` is a better solution
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

taxID2name <- function(taxID) {
    kegg_taxa <- readRDS(system.file("extdata/kegg_taxa.rds",
        package = "clusterProfiler"))
    kegg_taxa$kegg.name[kegg_taxa$kegg.taxa == taxID]
}


removeEmptyEntry.list <- function(x) {
    notNA.idx <- unlist(lapply(x, function(i) !is.null(i) && !all(is.na(i))))
    x[notNA.idx]
}


GSEA_internal <- DOSE:::GSEA_internal
enricher_internal <- DOSE:::enricher_internal

globalVariables(".")
