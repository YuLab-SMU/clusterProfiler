##' filter GO enriched result at specific level
##'
##' 
##' @title gofilter
##' @param x output from enrichGO or compareCluster
##' @param level GO level
##' @return updated object
##' @export
##' @author Guangchuang Yu
gofilter <- function(x, level=4) {
    ID <- NULL ## to satisfy codetools for using ID in subset function
    ont <- get_go_ontology(x)
    terms <- getGOLevel(ont, level)
    if (is(x, "enrichResult")) {
        x@result %<>% subset(ID %in% terms)
    } else {
        ## should be 'compareClusterResult', otherwise get_go_ontology already throw error
        x@compareClusterResult %<>% subset(ID %in% terms)
    }
    return(x)
}



