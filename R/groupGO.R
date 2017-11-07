
##' Functional Profile of a gene set at specific GO level.
##' Given a vector of genes, this function will return the GO profile at
##' a specific level.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param OrgDb OrgDb
##' @param keyType key type of input gene
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
groupGO <- function(gene, OrgDb, keyType="ENTREZID", ont="CC", level = 2, readable=FALSE) {
    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF"))

    GO_DATA <- get_GO_data(OrgDb, ont, keyType)

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
             keytype = keyType
             )
    if(readable == TRUE)
        x <- setReadable(x, OrgDb)

    return(x)
}

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
