#' simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms
#'
#'
#' @name simplify
#' @docType methods
#' @rdname simplify-methods
#' @title simplify method
#' @param x output of enrichGO
#' @param cutoff similarity cutoff
#' @param by feature to select representative term, selected by 'select_fun' function
#' @param select_fun function to select feature passed by 'by' parameter
#' @param measure method to measure similarity
#' @param semData GOSemSimDATA object
#' @return updated enrichResult object
#' @exportMethod simplify
#' @references issue #28
#' \url{https://github.com/GuangchuangYu/clusterProfiler/issues/28}
#' @aliases simplify,enrichResult-method
#' @author Guangchuang Yu
setMethod("simplify", signature(x="enrichResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData = NULL) {
              ontology <- x@ontology
              if (!ontology %in% c("BP", "MF", "CC", "GOALL")) {
                  ## `enricher()` (and other ORA entry points) may be run over a GO
                  ## gene-set collection without recording an ontology in @ontology.
                  ## Derive it from the GO IDs, so users no longer have to assign
                  ## the slot by hand (#369, #753).
                  ontology <- infer_go_ontology(as.data.frame(x)$ID)
              }
              if (!ontology %in% c("BP", "MF", "CC", "GOALL"))
                  stop("simplify only applied to output from gseGO and enrichGO, ",
                       "or to a result whose gene sets are GO terms...")
              res <- as.data.frame(x)
              if (ontology == "GOALL") {
                  x@result <- simplify_ALL(res = res, cutoff = cutoff, by = by,
                      select_fun = select_fun, measure = measure,
                      semData = semData)
              } else {
                  x@result <- simplify_internal(res = res, cutoff = cutoff,
                                    by = by, select_fun = select_fun, 
                                    measure = measure,
                                    ontology = ontology, 
                                    semData = semData)                      
              }
              return(x)
          }
)

#' @rdname simplify-methods
#' @exportMethod simplify
#' @references issue #162
#' \url{https://github.com/GuangchuangYu/clusterProfiler/issues/162}
#' @aliases simplify,gseaResult-method
#' @author Gwang-Jin Kim and Guangchuang Yu
setMethod("simplify", signature(x="gseaResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData=NULL) {
            ontology <- x@setType
            if (!ontology %in% c("BP", "MF", "CC", "GOALL")) {
                ## `simplify()` removes redundancy using GO semantic similarity, so
                ## it must know which ontology the terms belong to. A GSEA run over
                ## a GO gene-set collection (MSigDB C5, a gson GO file, ...) stores
                ## the collection name in @setType instead of an ontology, so derive
                ## it from the GO IDs themselves before giving up (#753).
                ontology <- infer_go_ontology(as.data.frame(x)$ID)
            }
            if (!ontology %in% c("BP", "MF", "CC", "GOALL")) {
                stop("simplify only applied to output from gseGO and enrichGO, ",
                     "or to a GSEA result whose gene sets are GO terms...")
            }
            res <- as.data.frame(x)
            if (ontology == "GOALL") {
                x@result <- simplify_ALL(res = res, cutoff = cutoff, by = by,
                                select_fun = select_fun, measure = measure,
                                semData = semData)
              } else {
                x@result <- simplify_internal(res = res, cutoff = cutoff,
                                by = by, select_fun = select_fun, 
                                measure = measure,
                                ontology = ontology, 
                                semData = semData)
              }
            return(x)
          }
)

#' Infer the GO ontology of a set of gene-set IDs
#'
#' Returns "BP", "CC" or "MF" when every ID is a GO term belonging to a single
#' ontology, "GOALL" when the GO terms span several ontologies, and the empty
#' string when the IDs are not GO terms (so the caller can keep its own error).
#'
#' @param ids character vector of gene-set IDs
#' @return a single string
#' @noRd
infer_go_ontology <- function(ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids)]
    if (!length(ids) || !all(grepl("^GO:[0-9]{7}$", ids))) {
        return("")
    }
    if (!requireNamespace("GO.db", quietly = TRUE)) {
        return("")
    }
    ont <- AnnotationDbi::Ontology(GO.db::GOTERM)[ids]
    ont <- unique(as.character(ont[!is.na(ont)]))
    if (length(ont) == 1) {
        return(ont)
    }
    if (length(ont) > 1) {
        return("GOALL")
    }
    ""
}

#' @importFrom GOSemSim mgoSim
#' @importFrom GOSemSim godata
#' @importFrom tidyr gather
simplify_internal <- function(res, cutoff=0.7, by="p.adjust", select_fun=min, 
                              measure="Rel", ontology, semData) {
    if (missing(semData) || is.null(semData)) {
        yulab.utils::yulab_warn("semData is not provided. It will be calculated automatically.")
        if (measure == "Wang") {
            semData <- godata(ont = ontology)
        } else {
            semData <- godata(ont = ontology, computeIC = TRUE)
        }
    } else {
        if (ontology != semData@ont) {
            msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
            stop(msg)
        }
    }

    sim <- mgoSim(res$ID, res$ID,
                  semData = semData,
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
    anc <- GOSemSim:::getAncestors(ontology)
    for (i in seq_along(ID)) {
        ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
        ## if length(ii) == 1, then go1 == go2
        if (length(ii) < 2)
            next

        sim_subset <- sim.df[ii,]

        jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))


        if (length(jj) > 1) {
            ll <- vapply(sim_subset$go1[jj], function(.id) length(anc[[.id]]), numeric(1))
            jj <- jj[which.max(ll)]
        }


        ## sim.df <- sim.df[-ii[-jj]]
        GO_to_remove <- unique(c(GO_to_remove, sim_subset$go1[-jj]))
    }

    res[!res$ID %in% GO_to_remove, ]
}



#' simplify output from compareCluster by removing redundancy of enriched GO terms
#'
#'
#' @name simplify
#' @docType methods
#' @rdname simplify-methods
#' @title simplify method
#' @return updated compareClusterResult object
#' @exportMethod simplify
#' @aliases simplify,compareClusterResult-method
#' @author Guangchuang Yu
setMethod("simplify", signature(x="compareClusterResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData=NULL) {
              res <- x@compareClusterResult
              ont <- get_go_ontology(x)

              ## organism <- x@.call$organism
              ## if (is.null(organism)) {
              ##     organism <- "human"
              ## }

              ## to satisfy codetools in calling subset
              Cluster <- NULL
              lres <- lapply(unique(res$Cluster), function(cls) subset(res, Cluster == cls))
              if (ont == "ALL") {
                  lres %<>% lapply(., simplify_ALL,
                      cutoff = cutoff, by = by,
                      select_fun = select_fun,
                      measure = measure,
                      semData = semData)
              } else {
                  lres %<>% lapply(., simplify_internal,
                                   cutoff=cutoff,
                                   by = by,
                                   select_fun = select_fun,
                                   measure = measure,
                                   ontology = ont,
                                   semData = semData)
              }
              x@compareClusterResult <- do.call("rbind", lres)
              return(x)
          }
)


#' @param data.frame of enrichment result 
#' @param cutoff similarity cutoff
#' @param by feature to select representative term, selected by 'select_fun' function
#' @param select_fun function to select feature passed by 'by' parameter
#' @param measure method to measure similarity
#' @param semData GOSemSimDATA object
#' @noRd
simplify_ALL <- function(res, cutoff, by, select_fun, measure, semData) {
    ONTOLOGY <- NULL
    lres <- lapply(unique(res[, "ONTOLOGY"]), function(y)
                      simplify_internal(dplyr::filter(res, ONTOLOGY == y),
                          cutoff = cutoff,
                          by = by,
                          select_fun = select_fun,
                          measure = measure,
                          ontology = y,
                          semData = NULL)
                  )
    do.call(rbind, lres)
}
