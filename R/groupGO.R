
##' Functional Profile of a gene set at specific GO level.
##' Given a vector of genes, this function will return the GO profile at
##' a specific level.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param OrgDb OrgDb
##' @param keytype keytype of input gene
##' @param ont One of "MF", "BP", and "CC" subontologies.
##' @param level Specific GO Level.
##' @param readable if readable is TRUE, the gene IDs will mapping to gene
##'   symbols.
##' @return A \code{groupGOResult} instance.
##' @seealso \code{\link{groupGOResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @importFrom methods new
##' @importClassesFrom methods data.frame
##' @importFrom DOSE setReadable
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @examples
##'
##' 	data(gcSample)
##' 	yy <- groupGO(gcSample[[1]], 'org.Hs.eg.db', ont="BP", level=2)
##' 	head(summary(yy))
##' 	#plot(yy)
##'
groupGO <- function(gene, OrgDb, keytype="ENTREZID", ont="CC", level = 2, readable=FALSE) {
    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF"))
    
    GO_DATA <- get_GO_data(OrgDb, ont, keytype)
    
    GOLevel <- getGOLevel(ont, level) ##get GO IDs of specific level.

    DOSE <- "DOSE"
    require(DOSE, character.only = TRUE)
    TERMID2EXTID <- eval(parse(text=paste0(DOSE, ":::", "TERMID2EXTID")))
    TERM2NAME <- eval(parse(text=paste0(DOSE, ":::", "TERM2NAME")))

    
    GO2ExtID <- TERMID2EXTID(GOLevel, GO_DATA) ## mapping GOID to External Gene IDs.

    geneID.list <- lapply(GO2ExtID, function(x) gene[gene %in% x]) ## retain External Gene IDs which appear in *gene*

    ## if (readable) {
        ## mapping Gene IDs to Gene Names.
    ##    geneID.list <- lapply(geneID.list, EXTID2NAME, organism=organism)
    ## }
    geneID <- sapply(geneID.list, function(i) paste(i, collapse="/"))

    Count <- unlist(lapply(geneID.list, length))
    GeneRatio <- paste(Count, length(unique(unlist(gene))), sep="/")
    Descriptions <- TERM2NAME(GOLevel, GO_DATA)
    result = data.frame(ID=as.character(GOLevel),
        Description=Descriptions,
        Count=Count,
        GeneRatio=GeneRatio,
        geneID=geneID)

    x <- new("groupGOResult",
             result=result,
             ontology = ont,
             level = level,
             organism = get_organism(OrgDb),
             gene = gene,
             keytype = keytype,
             geneInCategory = geneID.list
             )
    if(readable == TRUE) 
        x <- setReadable(x, OrgDb)
    
    return(x)
}

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
##' @slot geneInCategory gene and category association
##' @slot readable logical flag of gene ID in symbol or not.
##' @exportClass groupGOResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
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


## show method for \code{groupGOResult} instance
##
##
## @name show
## @docType methods
## @rdname show-methods
##
## @title show method
## @param object A \code{groupGOResult} instance
## @return message
## @importFrom methods show
## @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="groupGOResult"),
          function (object){
              ont = object@ontology
              Level = object@level
              Organism = object@organism
              Gene = object@gene
              cat ("GO", ont, "Profiles", "at level", Level, "of", length(Gene),  Organism, "genes", "\n")
          }
          )
