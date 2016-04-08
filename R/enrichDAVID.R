##' enrichment analysis by DAVID
##'
##' 
##' @title enrichDAVID
##' @param gene input gene
##' @param idType id type
##' @param listType list Type
##' @param minGSSize minimal size of genes annotated for testing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param annotation david annotation
##' @param pvalueCutoff pvalueCutoff
##' @param pAdjustMethod one of "BH" and "bonferroni"
##' @param qvalueCutoff qvalutCutoff
##' @param species species
##' @param david.user david user
##' @return A \code{enrichResult} instance
## @importFrom RDAVIDWebService DAVIDWebService
## @importFrom RDAVIDWebService addList
## @importFrom RDAVIDWebService setAnnotationCategories
## @importFrom RDAVIDWebService getFunctionalAnnotationChart
## @importFrom RDAVIDWebService getSpecieNames
##' @importFrom qvalue qvalue
##' @importFrom utils installed.packages
##' @export
##' @author Guangchuang Yu
enrichDAVID <- function(gene,
                        idType        = "ENTREZ_GENE_ID", 
                        listType      = "Gene",
                        minGSSize     = 10,
                        maxGSSize     = 500,
                        annotation    = "GOTERM_BP_FAT",
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.2,
                        species       = NA,
                        david.user){

    Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL
    
    pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", "BH"))
    idType <- match.arg(idType, c("AFFYMETRIX_3PRIME_IVT_ID",
                                  "AFFYMETRIX_EXON_GENE_ID",
                                  "AFFYMETRIX_SNP_ID",
                                  "AGILENT_CHIP_ID",
                                  "AGILENT_ID",
                                  "AGILENT_OLIGO_ID",
                                  "ENSEMBL_GENE_ID",
                                  "ENSEMBL_TRANSCRIPT_ID",
                                  "ENTREZ_GENE_ID",
                                  "GENOMIC_GI_ACCESSION",
                                  "GENPEPT_ACCESSION",
                                  "ILLUMINA_ID",
                                  "IPI_ID",
                                  "MGI_ID",
                                  "OFFICIAL_GENE_SYMBOL",
                                  "PFAM_ID",
                                  "PIR_ID",
                                  "PROTEIN_GI_ACCESSION",
                                  "REFSEQ_GENOMIC",
                                  "REFSEQ_MRNA",
                                  "REFSEQ_PROTEIN",
                                  "REFSEQ_RNA",
                                  "RGD_ID",
                                  "SGD_ID",
                                  "TAIR_ID",
                                  "UCSC_GENE_ID",
                                  "UNIGENE",
                                  "UNIPROT_ACCESSION",
                                  "UNIPROT_ID",
                                  "UNIREF100_ID",
                                  "WORMBASE_GENE_ID",
                                  "WORMPEP_ID",
                                  "ZFIN_ID"))

    david.pkg <- "RDAVIDWebService"
    pkgs <- installed.packages()[,1]
    if (! david.pkg %in% pkgs) {
        stop("You should have RDAVIDWebService package installed before using enrichDAVID...")
    }
    
    require(david.pkg, character.only=TRUE)
    DAVIDWebService <- eval(parse(text="DAVIDWebService"))
    addList <- eval(parse(text="addList"))
    setAnnotationCategories <- eval(parse(text="setAnnotationCategories"))
    getFunctionalAnnotationChart <- eval(parse(text="getFunctionalAnnotationChart"))
    getSpecieNames <- eval(parse(text="getSpecieNames"))

    
    david <- DAVIDWebService$new(email=david.user,
                                 url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    david.res <- addList(david, gene, idType=idType,
                         listName="clusterProfiler",
                         listType=listType)
                        
                        
    if (david.res$inDavid == 0) {
        stop("All id can not be mapped. Please check 'idType' parameter...")
    }
    
    setAnnotationCategories(david, annotation)
    x <- getFunctionalAnnotationChart(david, threshold=1, count=minGSSize)

    if (length(x@.Data) == 0) {
        warning("No significant enrichment found...")
        return(NULL)
    }
    
    term <- x$Term
    if (length(grep("~", term[1])) == 0) {
        sep <- ":"
    } else {
        sep <- "~"
    }
    term.list <- sapply(term, function(y) strsplit(y, split=sep))
    term.df <- do.call("rbind", term.list)
    ID <- term.df[,1]
    Description <- term.df[,2]
    GeneRatio <- with(x, paste(Count, List.Total, sep="/"))
    BgRatio <- with(x, paste(Pop.Hits, Pop.Total, sep="/"))
    Over <- data.frame(ID          = ID,
                       Description = Description,
                       GeneRatio   = GeneRatio,
                       BgRatio     = BgRatio,
                       pvalue      = x$PValue,
                       stringsAsFactors = FALSE)
    row.names(Over) <- ID

    if (pAdjustMethod == "bonferroni") {
        Over$p.adjust <- x$Bonferroni
    } else {
        Over$p.adjust <- x$Benjamini
    }

    qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"),
                     error=function(e) NULL)
    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    Over$qvalue <- qvalues
    Over$geneID <- gsub(",\\s*", "/", x$Genes)
    Over$Count <- x$Count

    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
    if (! any(is.na(Over$qvalue))) {
        Over <- Over[ Over$qvalue <= qvalueCutoff, ]
    }

    org <- getSpecieNames(david)
    org <- gsub("\\(.*\\)", "", org)

    gc <- strsplit(Over$geneID, "/")
    names(gc) <- Over$ID

    if (!is.na(maxGSSize) || !is.null(maxGSSize)) {
        idx <- x$Pop.Hits <= maxGSSize
        Over <- Over[idx,]
    }
    
    new("enrichResult",
        result         = Over,
        pvalueCutoff   = pvalueCutoff,
        pAdjustMethod  = pAdjustMethod,
        organism       = org,
        ontology       = annotation, ## as.character(x$Category[1]),
        gene           = as.character(gene),
        geneInCategory = gc,
        keytype        = idType)
}
    
