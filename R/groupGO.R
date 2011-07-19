

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
#' @examples
#'
#' 	data(gcSample)
#' 	yy <- groupGO(gcSample[[1]], organism="human", ont="BP", level=2)
#' 	head(summary(yy))
#' 	#plot(yy)
#'


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

setClass("groupGOResult",
         representation=representation(
         groupGOResult="data.frame",
         Ont = "character",
         Level = "numeric",
         Organism = "character",
         Gene = "character"
         )
         )

setMethod("show", signature(object="groupGOResult"),
          function (object){
              ont = object@Ont
              Level = object@Level
              Organism = object@Organism
              Gene = object@Gene
              cat ("GO", ont, "Profiles", "at level", Level, "of", length(Gene),  Organism, "genes", "\n")
          }
          )

setMethod("summary", signature(object="groupGOResult"),
          function (object){
              return(object@groupGOResult)
          }
          )

setMethod("plot", signature(x="groupGOResult"),
          function (x, order="FALSE", title="", font.size=12, drop=FALSE){
              groupGOResult <- summary(x)
              if (drop == TRUE) {
                                        # drop void category.
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
