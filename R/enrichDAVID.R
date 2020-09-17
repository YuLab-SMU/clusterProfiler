##' enrichment analysis by DAVID
##'
##'
##' @title enrichDAVID
##' @param gene input gene
##' @param idType id type
##' @param minGSSize minimal size of genes annotated for testing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param annotation david annotation
##' @inheritParams enricher
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
                        universe,
                        minGSSize     = 10,
                        maxGSSize     = 500,
                        annotation    = "GOTERM_BP_FAT",
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.2,
                        david.user){
    
    Count <- List.Total <- Pop.Hits <- Pop.Total <- NULL
    
    pAdjustMethod <- match.arg(pAdjustMethod, c("bonferroni", "BH"))
    
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
    getIdTypes <- eval(parse(text="getIdTypes"))
    
    david <- DAVIDWebService$new(email=david.user,
                                 url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    
    ## addList will throw error if idType is not match.
    ## use match.arg to check before addList make it more readable
    
    idType <- match.arg(idType, getIdTypes(david))
    
    ##     getIdTypes(david)
    ##  [1] "AFFYMETRIX_3PRIME_IVT_ID" "AFFYMETRIX_EXON_ID"      
    ##  [3] "AGILENT_CHIP_ID"          "AGILENT_ID"              
    ##  [5] "AGILENT_OLIGO_ID"         "APHIDBASE_ID"            
    ##  [7] "BEEBASE_ID"               "BEETLEBASE_ID"           
    ##  [9] "BGD_ID"                   "CGNC_ID"                 
    ## [11] "CRYPTODB_ID"              "DICTYBASE_ID"            
    ## [13] "ENSEMBL_GENE_ID"          "ENSEMBL_TRANSCRIPT_ID"   
    ## [15] "ENTREZ_GENE_ID"           "FLYBASE_GENE_ID"         
    ## [17] "GENBANK_ACCESSION"        "GENOMIC_GI_ACCESSION"    
    ## [19] "GENPEPT_ACCESSION"        "LOCUS_TAG"               
    ## [21] "MGI_ID"                   "MIRBASE_ID"              
    ## [23] "MRNA_GI_ACCESSION"        "NASONIABASE_ID"          
    ## [25] "PROTEIN_GI_ACCESSION"     "PSEUDOCAP_ID"            
    ## [27] "REFSEQ_MRNA"              "REFSEQ_PROTEIN"          
    ## [29] "RGD_ID"                   "SGD_ID"                  
    ## [31] "TAIR_ID"                  "UNIGENE"                 
    ## [33] "UNIPROT_ACCESSION"        "UNIPROT_ID"              
    ## [35] "VECTORBASE_ID"            "WORMBASE_GENE_ID"        
    ## [37] "XENBASE_ID"               "ZFIN_ID"
    
    david.res <- addList(david, gene, idType=idType,
                         listName="clusterProfiler",
                         listType="Gene")
    
    
    if (david.res$inDavid == 0) {
        stop("All id can not be mapped. Please check 'idType' parameter...")
    }
    
    if (!missing(universe)) {
        david.res <- addList(david, universe, idType=idType,
                             listName="universe",
                             listType="Background")
    }
    
    setAnnotationCategories(david, annotation)
    x <- getFunctionalAnnotationChart(david, threshold=1, count=minGSSize)
    
    if (!is.na(maxGSSize) && !is.null(maxGSSize)) {
        x <- x[x$Count <= maxGSSize, ]
    }
    
    if (nrow(x) == 0) {
        warning("No enriched GO term found...")
        return(NULL)
    }
    
    sep <- ifelse(grepl("~", x$Term[1]), "~", ":")
    term.df <- do.call("rbind", sapply(x$Term, function(y) strsplit(y, split=sep)))
    
    Over <- data.frame(ID          = term.df[,1],
                       Description = term.df[,2],
                       GeneRatio   = with(x, paste(Count, List.Total, sep="/")),
                       BgRatio     = with(x, paste(Pop.Hits, Pop.Total, sep="/")),
                       pvalue      = x$PValue,
                       p.adjust    = p.adjust(x$PValue, method = pAdjustMethod),
                       stringsAsFactors = FALSE)
    row.names(Over) <- Over$ID
    
    qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"),
                     error=function(e) NULL)
    if (class(qobj) == "qvalue") {
        Over$qvalue <- qobj$qvalues
    } else {
        Over$qvalue <- NA
    }
    
    Over$geneID <- gsub(",\\s*", "/", x$Genes)
    Over$Count <- x$Count
    
    Over <- Over[Over$p.adjust <= pvalueCutoff, ]
    
    if (! any(is.na(Over$qvalue))) {
        Over <- Over[Over$qvalue <= qvalueCutoff, ]
    }
    
    species <- getSpecieNames(david)
    species <- gsub("\\(.*\\)", "", species)
    
    new("enrichResult",
        result         = Over,
        pvalueCutoff   = pvalueCutoff,
        pAdjustMethod  = pAdjustMethod,
        organism       = org,
        ontology       = annotation,
        gene           = as.character(gene),
        keytype        = idType)
}

