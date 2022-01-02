##' Compare gene clusters functional profile
##'
##' Given a list of gene set, this function will compute profiles of each gene
##' cluster.
##'
##'
##' @param geneClusters a list of entrez gene id. Alternatively, a formula of type \code{Entrez~group}
##' or a formula of type \code{Entrez | logFC ~ group} for "gseGO", "gseKEGG" and "GSEA".
##' @param fun One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .
##' @param data if geneClusters is a formula, the data from which the clusters must be extracted.
##' @param ...  Other arguments.
##' @return A \code{clusterProfResult} instance.
##' @importFrom methods new
##' @importFrom stats formula
##' @importFrom plyr llply
##' @importFrom plyr ldply
##' @importFrom plyr dlply
##' @importFrom utils modifyList
##' @importClassesFrom DOSE compareClusterResult
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
##'   \code{\link{enrichGO}}
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
##' ## formula interface for GSEA algorithm
##' require(breastCancerMAINZ)
##' data(mainz)
##' require(stringr)
##' require("hgu133a.db")
##' require(plyr)
##' clmainz=pData(mainz)$grade
##' dd <- exprs(mainz)
##' g1 <- dd[,clmainz == 1]
##' g2 <- dd[,clmainz == 2]
##' g3 <- dd[,clmainz == 3]
##' 
##' buildGenelist <- function(mat1, mat2) {
##'     geneList <- exp(rowMeans(mat1))/exp(rowMeans(mat2))
##'     geneList <- sort(geneList, decreasing=TRUE)
##'     geneList <- log(geneList, base=2)
##'     
##'     eg <- mget(names(geneList), hgu133aENTREZID, ifnotfound=NA)   
##'     gg <- data.frame(probe=names(geneList), val = geneList)
##'     eg.df <- data.frame(probe=names(eg), eg=unlist(eg))
##'     
##'     xx <- merge(gg, eg.df, by.x="probe", by.y="probe")
##'     xx <- xx[,-1]
##'     xx <- unique(xx)
##'     xx <- xx[!is.na(xx[,2]),]
##'     yy <- ddply(xx, .(eg), function(x) data.frame(val=mean(x$val)))
##'     
##'     geneList <- yy$val
##'     names(geneList) <- yy$eg
##'     geneList <- sort(geneList, decreasing=TRUE)
##' }
##' 
##' geneList21 <- buildGenelist(g2, g1)
##' geneList31 <- buildGenelist(g3, g1)
##' geneList32 <- buildGenelist(g3, g2)
##' 
##' mydf2 <- data.frame(Entrez = c(names(geneList21), names(geneList31), names(geneList32)),
##'                     logFC = c(geneList21, geneList31, geneList32),
##'                     group = c(rep("grade2_1", length(geneList21)), 
##'                               rep("grade3_1", length(geneList31)), 
##'                               rep("grade3_2", length(geneList32))))
##' ## formula interface for gseGO
##' gsea.formula <- compareCluster(Entrez|logFC~group, data=mydf2,
##'                                fun='gseGO', OrgDb='org.Hs.eg.db',
##'                                eps = 0)
##' dotplot(gsea.formula)
##'                                  
##' ##  formula interface for gseKEGG                               
##' gsea.formula2 <- compareCluster(Entrez|logFC~group, data=mydf2,
##'                                 fun='gseKEGG', eps = 0)   
##' dotplot(gsea.formula2)
##' 
##' ## formula interface for GSEA
##' library(org.Hs.eg.db)
##' kegg_clu <- download_KEGG("hsa", keggType="KEGG")
##' TERM2GENE <- kegg_clu$KEGGPATHID2EXTID
##' TERM2NAME <- kegg_clu$KEGGPATHID2NAME
##'                                        
##' gsea.formula3 <- compareCluster(Entrez|logFC~group, data=mydf2,
##'                                 fun='GSEA', TERM2GENE = TERM2GENE, 
##'                                 TERM2NAME = TERM2NAME, eps = 0)    
##' dotplot(gsea.formula3)  
##' 
##' }
compareCluster <- function(geneClusters, fun="enrichGO", data='', ...) {

    if (is.character(fun)) {
        fun <- eval(parse(text=fun))
    }

    # Use formula interface for compareCluster
    if (typeof(geneClusters) == 'language') {
        if (!is.data.frame(data)) {
            stop ('no data provided with formula for compareCluster')
        } else {
            genes.var = all.vars(geneClusters)[1]
            n.var = length(all.vars(geneClusters))
            grouping.formula = gsub('^.*~', '~', as.character(as.expression(geneClusters)))   # For formulas like x~y+z
            n.group.var = length(all.vars(formula(grouping.formula)))
            geneClusters = dlply(.data=data, formula(grouping.formula), .fun=function(x) {
                if ( (n.var - n.group.var) == 1 ) {
                    as.character(x[[genes.var]])
                } else if ( (n.var - n.group.var) == 2 ) {
                    fc.var = all.vars(geneClusters)[2]
                    geneList = structure(x[[fc.var]], names = x[[genes.var]])
                    sort(geneList, decreasing=TRUE)
                } else {
                    stop('only Entrez~group or Entrez|logFC~group type formula is supported')
                }
            })
        }
    }
    clProf <- llply(geneClusters,
                    .fun=function(i) {
                        x=suppressMessages(fun(i, ...))
                        # if (class(x) == "enrichResult" || class(x) == "groupGOResult") {
                        if (inherits(x, c("enrichResult", "groupGOResult", "gseaResult"))){
                            as.data.frame(x)
                        }
                    }
                    )
    clusters.levels = names(geneClusters)
    clProf.df <- ldply(clProf, rbind)

    if (nrow(clProf.df) == 0) {
        warning("No enrichment found in any of gene cluster, please check your input...")
        return(NULL)
    }

    #clProf.df <- dplyr::rename(clProf.df, c(.id="Cluster"))
    clProf.df <- plyr::rename(clProf.df, c(.id="Cluster"))
    clProf.df$Cluster = factor(clProf.df$Cluster, levels=clusters.levels)

    if (is.data.frame(data) && grepl('+', grouping.formula)) {
        groupVarName <- strsplit(grouping.formula, split="\\+") %>% unlist %>%
            gsub("~", "", .) %>% gsub("^\\s*", "", .) %>% gsub("\\s*$", "", .)
        groupVars <- sapply(as.character(clProf.df$Cluster), strsplit, split="\\.") %>% do.call(rbind, .)
        for (i in seq_along(groupVarName)) {
            clProf.df[, groupVarName[i]] <- groupVars[,i]
        }
        i <- which(colnames(clProf.df) %in% groupVarName)
        j <- (1:ncol(clProf.df))[-c(1, i)]
        clProf.df <- clProf.df[, c(1, i, j)]
    }

    ##colnames(clProf.df)[1] <- "Cluster"
    res <- new("compareClusterResult",
               compareClusterResult = clProf.df,
               geneClusters = geneClusters,
               .call = match.call(expand.dots=TRUE)
               )

    params <- modifyList(extract_params(args(fun)),
                         extract_params(res@.call))

    keytype <- params[['keyType']]
    if (is.null(keytype)) keytype <- "UNKNOWN"
    readable <- params[['readable']]
    if (length(readable) == 0) readable <- FALSE
    
    res@keytype <- keytype
    res@readable <- as.logical(readable)
    res@fun <- params[['fun']]

    return(res)
}

extract_params <- function(x) {
    y <- rlang::quo_text(x)
    if (is.function(x)) y <- sub('\nNULL$', '', y)

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
## @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @importFrom utils str
setMethod("show", signature(object="compareClusterResult"),
          function (object){
              cmsg <- paste("T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, ",
                       "W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. ", 
                       "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ", 
                       "The Innovation. 2021, 2(3):100141",
                            sep="\n", collapse="\n")

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
              if (fun == "enrichDO" || fun == "enrichNCG") {
                  citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
                                        "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
                                        "  analysis. Bioinformatics 2015 31(4):608-609",
                                        sep="\n", collapse="\n")
              } else if (fun == "enrichPathway") {
                  citation_msg <- paste("  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
                                        "  reactome pathway analysis and visualization. Molecular BioSystems",
                                        "  2016, 12(2):477-479", sep="\n", collapse="\n")
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
## @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("summary", signature(object="compareClusterResult"),
          function(object, ...) {
              warning("summary method to convert the object to data.frame is deprecated, please use as.data.frame instead.")
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
##' @export
merge_result <- function(enrichResultList) {
    if ( !is(enrichResultList, "list")) {
        stop("input should be a name list...")
    }
    if ( is.null(names(enrichResultList))) {
        stop("input should be a name list...")
    }
    x <- lapply(enrichResultList, as.data.frame)
    names(x) <- names(enrichResultList)
    y <- ldply(x, "rbind")
    y <- plyr::rename(y, c(.id="Cluster"))
    y$Cluster = factor(y$Cluster, levels=names(enrichResultList))
    new("compareClusterResult",
        compareClusterResult = y)

}
