##' Class "compareClusterResult"
##' This class represents the comparison result of gene clusters by GO
##' categories at specific level or GO enrichment analysis.
##'
##'
##' @name compareClusterResult-class
##' @aliases compareClusterResult-class show,compareClusterResult-method
##'   summary,compareClusterResult-method plot,compareClusterResult-method
##' @docType class
##' @slot compareClusterResult cluster comparing result
##' @slot geneClusters a list of genes
##' @slot fun one of groupGO, enrichGO and enrichKEGG
##' @slot .call function call
##' @exportClass compareClusterResult
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @exportClass compareClusterResult
##' @seealso \code{\linkS4class{groupGOResult}}
##'   \code{\linkS4class{enrichResult}} \code{\link{compareCluster}}
##' @keywords classes
setClass("compareClusterResult",
         representation = representation(
             compareClusterResult = "data.frame",
             geneClusters = "list",
             fun = "character",
             .call = "call"
         )
         )


##' Class "groupGOResult"
##' This class represents the result of functional Profiles of a set of gene at
##' specific GO level.
##'
##'
##' @name groupGOResult-class
##' @aliases groupGOResult-class show,groupGOResult-method
##' @docType class
##' @slot result GO classification result
##' @slot ontology Ontology
##' @slot level GO level
##' @slot organism one of "human", "mouse" and "yeast"
##' @slot gene Gene IDs
##' @slot readable logical flag of gene ID in symbol or not.
##' @exportClass groupGOResult
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @importMethodsFrom DOSE summary
##' @importFrom DOSE setReadable
##' @seealso \code{\linkS4class{compareClusterResult}}
##'   \code{\link{compareCluster}} \code{\link{groupGO}}
##' @keywords classes
setClass("groupGOResult",
         representation=representation(
         level = "numeric"
         ),

         contains = "enrichResult"
         )
