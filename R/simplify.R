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
##' @importFrom GOSemSim mgoSim
##' @importFrom tidyr gather
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

              ## to satisfy codetools for calling gather
              go1 <- go2 <- similarity <- NULL
              
              res <- x@result
              
              sim <- mgoSim(res$ID, res$ID,
                            ont=x@ontology,
                            organism=x@organism,
                            measure=measure,
                            combine=NULL)
              
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
              
              x@result <- res[!res$ID %in% GO_to_remove, ]
              return(x)
          }
          )


