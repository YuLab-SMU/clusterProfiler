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
##' @author Guangchuang Yu \url{https://yulab-smu.top}
##' @importMethodsFrom DOSE summary
##' @importFrom DOSE setReadable
##' @seealso [compareClusterResult], [compareCluster], [groupGO]
##' @keywords classes
setClass(
  "groupGOResult",
  representation = representation(
    level = "numeric"
  ),

  contains = "enrichResult"
)
