##' Compare gene clusters functional profile
##'
##' Given a list of gene set, this function will compute profiles of each gene
##' cluster.
##'
##'
##' @param geneClusters a list of entrez gene id. Alternatively, a formula of type \code{Entrez~group}
##' or a formula of type \code{Entrez | logFC ~ group} for "gseGO", "gseKEGG" and "GSEA".
##' @param fun One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .
##' Users can also supply their own function.
##' @param data if geneClusters is a formula, the data from which the clusters must be extracted.
##' @param source_from If using a custom function in "fun", provide the source package as
##' a string here. Otherwise, the function will be obtained from the global environment.
##' @param ...  Other arguments.
##' @return A \code{clusterProfResult} instance.
##' @importFrom methods new
##' @importFrom stats formula
##' @importFrom plyr llply
##' @importFrom plyr ldply
##' @importFrom plyr dlply
##' @importFrom utils modifyList
##' @importFrom rlang %||%
##' @importClassesFrom DOSE compareClusterResult
##' @export
##' @author Guangchuang Yu \url{https://yulab-smu.top}
##' @seealso [compareClusterResult-class], [groupGO], [enrichGO], [enrichKEGG],
##' [enrichDO][DOSE::enrichDO], [enrichPathway][ReactomePA::enrichPathway]
##' @keywords manip
##' @examples
##' \dontrun{
##' data(gcSample)
##' xx <- compareCluster(gcSample, fun="enrichKEGG",
##'                      organism="hsa", pvalueCutoff=0.05)
##' as.data.frame(xx)
##' # plot(xx, type="dot", caption="KEGG Enrichment Comparison")
##' dotplot(xx)
##'
##' ## formula interface
##' mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
##'                             '100127206', '100128071'),
##'                    logFC = c(1.1, -0.5, 5, 2.5, -3, 3),
##'                    group = c('A', 'A', 'A', 'B', 'B', 'B'),
##'                    othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
##' xx.formula <- compareCluster(Entrez~group, data=mydf,
##'                              fun='groupGO', OrgDb='org.Hs.eg.db')
##' as.data.frame(xx.formula)
##'
##' ## formula interface with more than one grouping variable
##' xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf,
##'                                        fun='groupGO', OrgDb='org.Hs.eg.db')
##' as.data.frame(xx.formula.twogroups)
##'
##' }
compareCluster <- function(
    geneClusters,
    fun = "enrichGO",
    data = '',
    source_from = NULL,
    ...
) {
    if (is.character(fun)) {
        if (
            fun %in%
                c(
                    "groupGO",
                    "enrichGO",
                    "enrichKEGG",
                    "gseGO",
                    "gseKEGG",
                    "GSEA",
                    "gseWP"
                )
        ) {
            fun <- utils::getFromNamespace(fun, "clusterProfiler")
        } else if (
            fun %in%
                c(
                    "enrichDO",
                    "enrichDGN",
                    "enrichDGNv",
                    "enrichNCG",
                    "gseDO",
                    "gseNCG",
                    "gseDGN"
                )
        ) {
            check_installed(
                'DOSE',
                paste0('for compareCluster with "fun =', fun, ' ."'),
                action = BiocManager::install
            )
            fun <- utils::getFromNamespace(fun, "DOSE")
        } else if (fun %in% c("enrichPathway", "gsePathway")) {
            check_installed(
                'ReactomePA',
                paste0('for compareCluster with "fun =', fun, ' ."'),
                action = BiocManager::install
            )
            fun <- utils::getFromNamespace(fun, "ReactomePA")
        } else if (fun %in% c("enrichMeSH", "gseMeSH")) {
            check_installed(
                'meshes',
                paste0('for compareCluster with "fun =', fun, ' ."'),
                action = BiocManager::install
            )
            fun <- utils::getFromNamespace(fun, "meshes")
        } else {
            source_env <- .GlobalEnv
            if (!is.null(source_from)) {
                source_env <- loadNamespace(source_from)
            }
            # If fun is in global or any loaded package, this will get it
            # This assumes that a user will actually load said package.
            fun <- get(fun, envir = source_env)
        }
    }

    # Use formula interface for compareCluster
    if (typeof(geneClusters) == 'language') {
        if (!is.data.frame(data)) {
            stop('no data provided with formula for compareCluster')
        } else {
            genes.var = all.vars(geneClusters)[1]
            n.var = length(all.vars(geneClusters))
            # For formulas like x~y+z
            grouping.formula = gsub(
                '^.*~',
                '~',
                as.character(as.expression(geneClusters))
            )
            n.group.var = length(all.vars(formula(grouping.formula)))
            geneClusters = dlply(
                .data = data,
                formula(grouping.formula),
                .fun = function(x) {
                    if ((n.var - n.group.var) == 1) {
                        as.character(x[[genes.var]])
                    } else if ((n.var - n.group.var) == 2) {
                        fc.var = all.vars(geneClusters)[2]
                        geneList = structure(
                            x[[fc.var]],
                            names = x[[genes.var]]
                        )
                        sort(geneList, decreasing = TRUE)
                    } else {
                        stop(
                            'only Entrez~group or Entrez|logFC~group type formula is supported'
                        )
                    }
                }
            )
        }
    }
    clProf <- llply(geneClusters, .fun = function(i) {
        x = suppressMessages(fun(i, ...))

        if (inherits(x, c("enrichResult", "groupGOResult", "gseaResult"))) {
            as.data.frame(x)
        }
    })
    clusters.levels = names(geneClusters)
    clProf.df <- ldply(clProf, rbind)

    if (nrow(clProf.df) == 0) {
        warning(
            "No enrichment found in any of gene cluster, please check your input..."
        )
        return(NULL)
    }

    #clProf.df <- dplyr::rename(clProf.df, c(.id="Cluster"))
    clProf.df <- plyr::rename(clProf.df, c(.id = "Cluster"))
    clProf.df$Cluster = factor(clProf.df$Cluster, levels = clusters.levels)

    if (is.data.frame(data) && grepl('+', grouping.formula)) {
        groupVarName <- strsplit(grouping.formula, split = "\\+") %>%
            unlist %>%
            gsub("~", "", .) %>%
            gsub("^\\s*", "", .) %>%
            gsub("\\s*$", "", .)
        groupVars <- sapply(
            as.character(clProf.df$Cluster),
            strsplit,
            split = "\\."
        ) %>%
            do.call(rbind, .)
        for (i in seq_along(groupVarName)) {
            clProf.df[, groupVarName[i]] <- groupVars[, i]
        }
        i <- which(colnames(clProf.df) %in% groupVarName)
        j <- (1:ncol(clProf.df))[-c(1, i)]
        clProf.df <- clProf.df[, c(1, i, j)]
    }

    ##colnames(clProf.df)[1] <- "Cluster"
    res <- new(
        "compareClusterResult",
        compareClusterResult = clProf.df,
        geneClusters = geneClusters,
        .call = match.call(expand.dots = TRUE)
    )

    params <- modifyList(extract_params(args(fun)), extract_params(res@.call))

    keytype <- params[['keyType']]
    if (is.null(keytype)) {
        keytype <- "UNKNOWN"
    }
    readable <- params[['readable']]
    if (length(readable) == 0) {
        readable <- FALSE
    }

    res@keytype <- keytype
    res@readable <- as.logical(readable)
    ## work-around for bug in extract_parameters -- it doesn't match default args
    res@fun <- params[['fun']] %||% 'enrichGO'

    return(res)
}

extract_params <- function(x) {
    y <- rlang::quo_text(x)
    if (is.function(x)) {
        y <- sub('\nNULL$', '', y)
    }

    y <- gsub('"', '', y) %>%
        ## sub(".*\\(", "", .) %>%
        sub("[^\\(]+\\(", "", .) %>%
        sub("\\)$", "", .) %>%
        gsub("\\s+", "", .)

    y <- strsplit(y, ",")[[1]]
    params <- sub("=.*", "", y)
    vals <- sub(".*=", "", y)
    i <- params != vals
    params <- params[i]
    vals <- vals[i]
    names(vals) <- params
    return(as.list(vals))
}


## show method for \code{compareClusterResult} instance
##
##
## @name show
## @alias show
## @docType methods
## @rdname show-methods
##
## @title show method
## @param object A \code{compareClusterResult} instance.
## @return message
## @importFrom methods show
## @author Guangchuang Yu \url{https://yulab-smu.top}
##' @importFrom utils str
setMethod("show", signature(object = "compareClusterResult"), function(object) {
    cmsg <- paste(
        "T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, ",
        "W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. ",
        "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ",
        "The Innovation. 2021, 2(3):100141",
        sep = "\n",
        collapse = "\n"
    )

    geneClusterLen <- length(object@geneClusters)
    fun <- object@fun
    result <- object@compareClusterResult
    clusts <- split(result, result$Cluster)
    nterms <- sapply(clusts, nrow)

    cat("#\n# Result of Comparing", geneClusterLen, "gene clusters", "\n#\n")
    cat("#.. @fun", "\t", fun, "\n")
    cat("#.. @geneClusters", "\t")
    str(object@geneClusters)
    cat("#...Result", "\t")
    str(result)
    cat("#.. number of enriched terms found for each gene cluster:\n")
    for (i in seq_along(clusts)) {
        cat("#..  ", paste0(names(nterms)[i], ":"), nterms[i], "\n")
    }
    cat("#\n#...Citation\n")
    citation_msg <- NULL
    if (length(fun) == 0) {
        # do nothing
    } else if (fun == "enrichDO" || fun == "enrichNCG") {
        citation_msg <- paste(
            "  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
            "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
            "  analysis. Bioinformatics 2015 31(4):608-609",
            sep = "\n",
            collapse = "\n"
        )
    } else if (fun == "enrichPathway") {
        citation_msg <- paste(
            "  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
            "  reactome pathway analysis and visualization. Molecular BioSystems",
            "  2016, 12(2):477-479",
            sep = "\n",
            collapse = "\n"
        )
    }
    if (!is.null(citation_msg)) {
        cat(paste0("1.", citation_msg), "\n\n")
        cat(paste0("2.", cmsg), "\n\n")
    } else {
        cat(cmsg, "\n\n")
    }
})

## summary method for \code{compareClusterResult} instance
##
##
## @name summary
## @alias summary
## @docType methods
## @rdname summary-methods
##
## @title summary method
## @param object A \code{compareClusterResult} instance.
## @return A data frame
## @importFrom stats4 summary
## @exportMethod summary
## @author Guangchuang Yu \url{https://yulab-smu.top}
setMethod(
    "summary",
    signature(object = "compareClusterResult"),
    function(object, ...) {
        warning(
            "summary method to convert the object to data.frame is deprecated, please use as.data.frame instead."
        )
        return(as.data.frame(object, ...))
    }
)


##' merge a list of enrichResult objects to compareClusterResult
##'
##'
##' @title merge_result
##' @param enrichResultList a list of enrichResult objects
##' @return a compareClusterResult instance
##' @author Guangchuang Yu
##' @importFrom plyr ldply
##' @importFrom methods is
##' @export
merge_result <- function(enrichResultList) {
    if (!is(enrichResultList, "list")) {
        stop("input should be a name list...")
    }
    if (is.null(names(enrichResultList))) {
        stop("input should be a name list...")
    }
    x <- lapply(enrichResultList, as.data.frame)
    names(x) <- names(enrichResultList)
    y <- ldply(x, "rbind")
    y <- plyr::rename(y, c(.id = "Cluster"))
    y$Cluster = factor(y$Cluster, levels = names(enrichResultList))
    new("compareClusterResult", compareClusterResult = y)
}
