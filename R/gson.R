##' download the latest version of KEGG pathway and stored in a 'GSON' object
##'
##'
##' @title gson_KEGG
##' @param species species
##' @return a 'GSON' object
##' @author Guangchuang Yu
##' @importFrom gson gson
##' @export
gson_KEGG <- function(species) {
    x <- download_KEGG(species)
    gsid2gene <- setNames(x[[1]], c("gsid", "gene"))
    gsid2name <- setNames(x[[2]], c("gsid", "name"))
    y <- readLines("https://rest.kegg.jp/info/hsa")
    version <- sub("\\w+\\s+", "", y[grep('Release', y)])
    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name,
        species = species,
        version = version,
        accessed_date = as.character(Sys.Date())
    )
}

