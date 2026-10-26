#' simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms
#'
#' `simplify()` removes redundant enriched terms by grouping the terms whose
#' similarity exceeds `cutoff` and keeping one representative per group.
#'
#' For GO results the similarity is semantic and comes from GOSemSim
#' (`measure = "Wang"`, `"Resnik"`, `"Rel"`, `"Jiang"`, `"Lin"`, `"TC"`),
#' which needs an ontology and a `GOSemSimDATA` object.
#'
#' For KEGG, Reactome, MSigDB and any other gene-set collection, set
#' `measure = "JC"` (Jaccard index of the gene sets) or
#' `measure = "overlap"` (shared genes over the smaller set). These measures
#' are computed from `@geneSets`, need no ontology, and therefore work on any
#' enrichment result, including Reactome and KEGG results for which GO
#' semantic similarity is undefined (#743). A Reactome or KEGG hierarchy gives
#' many significantly enriched terms that share most of their genes; keeping
#' the best-scoring term of each overlapping group compresses them into the
#' distinct biological programs.
#'
#' @name simplify
#' @docType methods
#' @rdname simplify-methods
#' @title simplify method
#' @param x output of enrichGO
#' @param cutoff similarity cutoff
#' @param by feature to select representative term, selected by 'select_fun' function
#' @param select_fun function to select feature passed by 'by' parameter
#' @param measure method to measure similarity. One of `"Wang"` (default),
#'   `"Resnik"`, `"Rel"`, `"Jiang"`, `"Lin"`, `"TC"` for GO semantic
#'   similarity, or `"JC"` / `"overlap"` to compare the term gene sets, which
#'   works for every collection (KEGG, Reactome, MSigDB, ...)
#' @param semData GOSemSimDATA object
#' @return updated enrichResult object
#' @exportMethod simplify
#' @references issue #28
#' \url{https://github.com/GuangchuangYu/clusterProfiler/issues/28}
#' @seealso [enrichplot::pairwise_termsim]
#' @aliases simplify,enrichResult-method
#' @author Guangchuang Yu
setMethod("simplify", signature(x="enrichResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData = NULL) {
              if (is_gene_overlap_measure(measure)) {
                  x@result <- simplify_gene_overlap(x,
                      cutoff = cutoff, by = by, select_fun = select_fun,
                      measure = measure)
                  return(x)
              }

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
                       "or to a result whose gene sets are GO terms... ",
                       "For KEGG, Reactome and other non-GO results, set ",
                       "measure = \"JC\" to remove redundancy by gene overlap.")
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
            if (is_gene_overlap_measure(measure)) {
                x@result <- simplify_gene_overlap(x,
                    cutoff = cutoff, by = by, select_fun = select_fun,
                    measure = measure)
                return(x)
            }

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
                     "or to a GSEA result whose gene sets are GO terms... ",
                     "For KEGG, Reactome and other non-GO results, set ",
                     "measure = \"JC\" to remove redundancy by gene overlap.")
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


## Gene-overlap based redundancy removal -------------------------------
##
## `simplify()` historically measured redundancy with GO semantic
## similarity, which needs an ontology and a GOSemSimDATA object. KEGG,
## Reactome, MSigDB and every other non-GO collection returned results that
## `simplify()` refused to touch (#743). Those collections are already
## gene-set annotations, so the genes themselves give the similarity:
## measure = "JC" compares the term's gene set with the gene sets of the
## other terms, exactly like the emapplot similarity used by enrichplot.

#' measures that are computed from gene sets rather than a GO ontology
#' @noRd
is_gene_overlap_measure <- function(measure) {
    is.character(measure) && length(measure) == 1L &&
        !is.na(measure) &&
        toupper(measure) %in% c("JC", "JACCARD", "OVERLAP")
}

#' gene sets of the enriched terms, for overlap-based similarity
#'
#' Terms whose gene set is missing from `@geneSets` are reconstructed from the
#' `geneID` / `core_enrichment` column of the result, so results that were read
#' back from a table still work.
#' @noRd
gene_overlap_sets <- function(x, res) {
    ids <- as.character(res$ID)
    gene_sets <- x@geneSets

    sets <- lapply(ids, function(id) {
        gs <- gene_sets[[id]]
        if (!is.null(gs)) {
            return(unique(as.character(gs)))
        }
        col <- intersect(c("geneID", "core_enrichment"), names(res))
        if (length(col) == 0L) {
            return(character(0))
        }
        unique(strsplit(as.character(res[[col[1]]][match(id, ids)]),
                        "/", fixed = TRUE)[[1]])
    })
    names(sets) <- ids
    sets
}

#' Jaccard (or overlap) similarity matrix of gene sets
#'
#' Both measures are symmetric, so only the upper triangle is computed.
#' @noRd
gene_overlap_similarity <- function(sets, measure) {
    n <- length(sets)
    sim <- matrix(0, n, n)
    if (n == 0L) {
        return(sim)
    }
    if (n == 1L) {
        sim[1, 1] <- 1
        return(sim)
    }
    jaccard <- !toupper(measure) %in% c("OVERLAP")
    for (i in seq_len(n - 1L)) {
        a <- sets[[i]]
        if (length(a) == 0L) {
            next
        }
        for (j in seq.int(i + 1L, n)) {
            b <- sets[[j]]
            if (length(b) == 0L) {
                next
            }
            shared <- length(intersect(a, b))
            if (shared == 0L) {
                next
            }
            s <- if (jaccard) {
                shared / length(union(a, b))
            } else {
                shared / min(length(a), length(b))
            }
            sim[i, j] <- s
            sim[j, i] <- s
        }
    }
    diag(sim) <- 1
    sim
}

#' remove redundant terms by gene-set overlap
#'
#' Same greedy logic as `simplify_internal()`: for every term, the terms whose
#' similarity exceeds `cutoff` form a redundant group, and only the
#' representative picked by `select_fun(by)` is kept. Ties are broken towards
#' the larger gene set when `select_fun` is `min`, the way the GO version
#' breaks ties towards the more specific term.
#' @noRd
simplify_gene_overlap <- function(x, cutoff = 0.7, by = "p.adjust",
                                  select_fun = min, measure = "JC") {
    if (!is_gene_overlap_measure(measure)) {
        measure <- "JC"
    }
    res <- as.data.frame(x)
    if (nrow(res) == 0L) {
        return(res)
    }
    if (!by %in% colnames(res)) {
        stop("'by' should be one of the columns of the result: ",
             paste(colnames(res), collapse = ", "), call. = FALSE)
    }

    ids <- as.character(res$ID)
    sets <- gene_overlap_sets(x, res)
    sim <- gene_overlap_similarity(sets, measure)
    dimnames(sim) <- list(ids, ids)

    by_value <- res[[by]]
    set_size <- vapply(sets, length, integer(1))
    ## prefer the bigger gene set on ties, and the earlier row as a last resort
    tie_key <- -set_size
    if (identical(select_fun, min)) {
        rank_secondary <- order(by_value, tie_key, seq_along(ids), na.last = TRUE)
    } else {
        rank_secondary <- order(by_value, seq_along(ids), na.last = TRUE)
    }
    priority <- match(seq_along(ids), rank_secondary)

    removed <- logical(nrow(res))
    for (k in seq_along(ids)) {
        if (removed[k]) {
            next
        }
        group <- which(sim[k, ] > cutoff)
        if (length(group) < 2L) {
            next
        }
        winner <- group[which.min(priority[group])]
        removed[group[group != winner]] <- TRUE
    }

    res[!removed, , drop = FALSE]
}
