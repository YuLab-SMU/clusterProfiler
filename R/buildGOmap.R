##' building GO mapping files
##'
##' provided by a data.frame of GO (column 1) and gene (column 2) direct annotation
##' this function will building gene to GO and GO to gene mapping,
##' with directly and undirectly (ancestor GO term) annotation.
##' @title buildGOmap
##' @param gomap data.frame with two columns of GO and gene ID
##' @return data.frame, GO annotation with indirect annotation
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
##' @importFrom utils stack
##' @export
##' @author Yu Guangchuang
buildGOmap <- function(gomap) {

    ## remove empty GO annotation
    gomap <- gomap[gomap[,1] != "", ]
    
    Gene2GO <- split(as.character(gomap[,1]), as.character(gomap[,2]))
    
    Gene2ALLGO <- lapply(Gene2GO,
                         function(i) {
                             mfans <- unlist(mget(i, GOMFANCESTOR, ifnotfound=NA))
                             bpans <- unlist(mget(i, GOBPANCESTOR, ifnotfound=NA))
                             ccans <- unlist(mget(i, GOCCANCESTOR, ifnotfound=NA))
                             ans <- c(mfans, bpans, ccans)
                             ans <- ans[ !is.na(ans) ]
                             ans <- c(i, ans)
                             ans <- unique(ans)
                             ans <- ans[ans != "all"]
                             return(ans)
                         })

    go2gene <- stack(Gene2ALLGO)
    colnames(go2gene) <- c("GO", "Gene")
    
    return(go2gene)
}
