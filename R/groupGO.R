#' Functional Profile of a gene set at specific GO level.
#' Given a vector of genes, this function will return the GO profile at
#' specific level.
#'
#'
#' @param gene a vector of entrez gene id.
#' @param organism Currently, only "human" and "mouse" supported.
#' @param ont One of "MF", "BP", and "CC" subontologies.
#' @param level Specific GO Level.
#' @param readable if readable is TRUE, the gene IDs will mapping to gene
#'   symbols.
#' @return A \code{groupGOResult} instance.
#' @seealso \code{\link{groupGOResult-class}}, \code{\link{compareCluster}}
#' @keywords manip
#' @export
#' @examples
#'
#' 	data(gcSample)
#' 	yy <- groupGO(gcSample[[1]], organism="human", ont="BP", level=2)
#' 	head(summary(yy))
#' 	#plot(yy)
#'
groupGO <- function(gene, organism="human", ont="CC", level = 2, readable=FALSE) {
    GOLevel <- getGOLevel(ont, level) ##get GO IDs of specific level.

    GO2ExtID <- getGO2ExtID(GOLevel, organism) ## mapping GOID to External Gene IDs.

    geneID.list <- lapply(GO2ExtID, function(x) gene[gene %in% x]) ## retain External Gene IDs which appear in *gene*

    if (readable) {
        geneID.list <- geneID2geneName(geneID.list, organism) ## mapping Gene IDs to Gene Names.
    }
    geneID <- sapply(geneID.list, function(i) paste(i, collapse="/"))

    Count <- unlist(lapply(geneID.list, length))
    Descriptions <- GO2Term(GOLevel)
    result = data.frame(GOID=GOLevel, Description=Descriptions, Count=Count, GeneID=geneID)
    new("groupGOResult",
        groupGOResult=result,
        Ont = ont,
        Level = level,
        Organism = organism,
        Gene = gene
	)
}

##' An S4 class that stores Gene Ontology classification result
##' @slot groupGOResult GO classification result
##' @slot Ont Ontology
##' @slot Level GO level
##' @slot Organism one of "human", "mouse" and "yeast"
##' @slot Gene Gene IDs
##' @author Guangchuang Yu
setClass("groupGOResult",
         representation=representation(
         groupGOResult="data.frame",
         Ont = "character",
         Level = "numeric",
         Organism = "character",
         Gene = "character"
         )
         )

##' show method for \code{groupGOResult} instance
##'
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{groupGOResult} instance
##' @return message
##' @author Guangchuang Yu
setMethod("show", signature(object="groupGOResult"),
          function (object){
              ont = object@Ont
              Level = object@Level
              Organism = object@Organism
              Gene = object@Gene
              cat ("GO", ont, "Profiles", "at level", Level, "of", length(Gene),  Organism, "genes", "\n")
          }
          )

##' summary method for \code{groupGOResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{groupGOResult} instance
##' @return A data frame
##' @author Guangchuang Yu
setMethod("summary", signature(object="groupGOResult"),
          function (object){
              return(object@groupGOResult)
          }
          )

##' plot method for \code{groupGOResult} instance
##'
##'
##' @name plot
##' @docType methods
##' @rdname plot-methods
##'
##' @title plot method
##' @param x A \code{groupGOResult} instance
##' @param order logical parameter, order the result by *Count*.
##' @param title graph title
##' @param font.size graph font size
##' @param drop logical parameter, drop void category.
##' @return ggplot object
##' @author GuangchuangYu
setMethod("plot", signature(x="groupGOResult"),
          function (x, order="FALSE", title="", font.size=12, drop=FALSE){
              groupGOResult <- summary(x)
              if (drop == TRUE) {
                  groupGOResult <- groupGOResult[groupGOResult$Count != 0, ]
              }
              if (order == TRUE) {
                  idx <- order(groupGOResult$Count)
                  groupGOResult <- groupGOResult[idx,]
              }
              groupGOResult$Description <- factor(groupGOResult$Description, level= as.character(groupGOResult$Description))
              p <- .barplotInternal(groupGOResult, title, font.size)
              p <- p +
                  aes(fill=Description) +
                      opts(legend.position="none")
              print(p)
          }
          )
