##' simplify output from enrichGO by removing redundancy of enriched GO terms
##'
##'
##' @name simplify
##' @docType methods
##' @rdname simplify-methods
##' @title simplify method
##' @param x output of enrichGO
##' @param cutoff similarity cutoff
##' @param by feature to select representative term, selected by 'select_fun' function
##' @param select_fun function to select feature passed by 'by' parameter
##' @param measure method to measure similarity
##' @importFrom IRanges simplify
##' @return updated enrichResult object
##' @exportMethod simplify
##' @references issue #28
##' \url{https://github.com/GuangchuangYu/clusterProfiler/issues/28}
##' @aliases simplify,enrichResult-method
##' @author Guangchuang Yu
setMethod("simplify", signature(x="enrichResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel") {
              if (!x@ontology %in% c("BP", "MF", "CC"))
                  stop("simplify only applied to output from enrichGO...")


              x@result %<>% simplify_internal(., cutoff, by, select_fun,
                                              measure, x@ontology, x@organism)

              return(x)
          }
          )

##' @importFrom GOSemSim mgoSim
##' @importFrom tidyr gather
simplify_internal <- function(res, cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel", ont, organism) {
    sim <- mgoSim(res$ID, res$ID,
                  ont=ont, 
                  organism=organism,
                  measure=measure,
                  combine=NULL)

    ## to satisfy codetools for calling gather
    go1 <- go2 <- similarity <- NULL
    
    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- gather(sim.df, go2, similarity, -go1)
    
    sim.df <- sim.df[!is.na(sim.df$similarity),]
    
    ## feature 'by' is attached to 'go1'
    sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
    sim.df$go2 <- as.character(sim.df$go2)
    
    ID <- res$ID
    
    GO_to_remove <- character()
    for (i in seq_along(ID)) {
        ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
        ## if length(ii) == 1, then go1 == go2
        if (length(ii) < 2) 
            next
        
        sim_subset <- sim.df[ii,]
        
        jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
        
        ## sim.df <- sim.df[-ii[-jj]]
        GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
    }
    
    res[!res$ID %in% GO_to_remove, ]
}


##' simplify output from compareCluster by removing redundancy of enriched GO terms
##'
##'
##' @name simplify
##' @docType methods
##' @rdname simplify-methods
##' @title simplify method
##' @importFrom IRanges simplify
##' @return updated compareClusterResult object
##' @exportMethod simplify
##' @aliases simplify,compareClusterResult-method
##' @author Guangchuang Yu
setMethod("simplify", signature(x="compareClusterResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel") {
              res <- x@compareClusterResult
              if (x@fun != "enrichGO") {
                  stop("simplify only work for GO...")
              }

              ont <- x@.call$ont
              if (is.null(ont)) {
                  ## should be "MF", default value of enrichGO
                  ## it's safe to determine from the output
                  ont <- res$ID[1] %>% GOTERM[[.]] %>% Ontology
              }

              organism <- x@.call$organism
              if (is.null(organism)) {
                  organism <- "human"
              }

              ## to satisfy codetools in calling subset
              Cluster <- NULL
              lres <- lapply(unique(res$Cluster), function(cls) subset(res, Cluster == cls)) 

              lres %<>% lapply(., simplify_internal,
                               cutoff=cutoff,
                               by = by,
                               select_fun = select_fun,
                               measure = measure,
                               ont = ont,
                               organism = organism)
              
              x@compareClusterResult <- do.call("rbind", lres)
              return(x)
          }
          )


