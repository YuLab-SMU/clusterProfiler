##' ORA analysis for WikiPathways
##'
##' This function performs over-representation analysis using WikiPathways
##' @title enrichWP
##' @param gene a vector of entrez gene id
##' @param organism a supported organism, which can be accessed via 
##' the get_wp_organisms() function, or a GSON object.
##' @param ... additional parameters, see also the parameters supported by the enricher() function
##' @return A \code{enrichResult} instance
##' @export
##' @author Guangchuang Yu
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' gene <- names(geneList)[abs(geneList) > 2]
##' ewp <- enrichWP(gene, organism = "Homo sapiens") 

##' wp_gson <- gson_WP("Homo sapiens")
##' ewp2 <- enrichWP(gene, organism = wp_gson) 
##' }
enrichWP <- function(gene, organism, ...) {


    if (inherits(organism, "character")) {                       
        wpdata <- prepare_WP_data(organism)
        TERM2GENE = wpdata$WPID2GENE
        TERM2NAME = wpdata$WPID2NAME
        species <- organism
    } else if (inherits(organism, "GSON")) {
        TERM2GENE = organism@gsid2gene
        TERM2NAME = organism@gsid2name
        species <- organism@species
    } else {
        stop("organism should be a species name or a GSON object")
    }

    res <- enricher(gene,
                    TERM2GENE = TERM2GENE,
                    TERM2NAME = TERM2NAME,
                    ...)
    if (is.null(res)) return(res)

    res@ontology <- "WikiPathways"
    res@organism <- species
    res@keytype <-  "ENTREZID"

    return(res)
}


##' GSEA analysis for WikiPathways
##'
##' This function performs GSEA using WikiPathways
##' @title gseWP
##' @param geneList ranked gene list
##' @param organism a supported organism, which can be accessed via the get_wp_organisms() function, or a GSON object.
##' @param ... additional parameters, see also the parameters supported by the GSEA() function
##' @return A \code{gseaResult} instance
##' @export
##' @author Guangchuang Yu 
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' gsewp <- gseWP(geneList, organism = "Homo sapiens") 

##' wp_gson <- gson_WP("Homo sapiens")
##' gsewp2 <- gseWP(geneList, organism = wp_gson) 
##' }
gseWP <- function(geneList, organism, ...) {

    
    if (inherits(organism, "character")) {                       
        wpdata <- prepare_WP_data(organism)
        TERM2GENE = wpdata$WPID2GENE
        TERM2NAME = wpdata$WPID2NAME
        species <- organism
    } else if (inherits(organism, "GSON")) {
        TERM2GENE = organism@gsid2gene
        TERM2NAME = organism@gsid2name
        species <- organism@species
    } else {
        stop("organism should be a species name or a GSON object")
    }

    res <- GSEA(geneList,
                TERM2GENE = TERM2GENE,
                TERM2NAME = TERM2NAME,
                ...)

    if (is.null(res)) return(res)

    res@setType <- "WikiPathways"
    res@organism <- species
    res@keytype <-  "ENTREZID"

    return(res)
}

##' @importFrom rlang .data
prepare_WP_data <- function(organism) {
    wp2gene <- get_wp_data(organism)
    ##TERM2GENE
    wpid2gene <- wp2gene %>% dplyr::select(.data$wpid, .data$gene) 
    ##TERM2NAME
    wpid2name <- wp2gene %>% dplyr::select(.data$wpid, .data$name) 
    list(WPID2GENE = wpid2gene,
         WPID2NAME = wpid2name)
}

get_wp_gmtfile <- function() {
    wpurl <- 'https://data.wikipathways.org/current/gmt/'
    x <- readLines(wpurl)
    y <- x[grep('\\.gmt',x)]
    sub(".*(wikipathways-.*\\.gmt).*", "\\1",  y[grep('File', y)])
}


##' list supported organism of WikiPathways
##'
##' This function extracts information from 'https://data.wikipathways.org/current/gmt/'
##' and lists all supported organisms
##' @title get_wp_organism
##' @return supported organism list
##' @export
##' @author Guangchuang Yu
get_wp_organisms <- function() {
    gmtfile <- get_wp_gmtfile()
    orgs <- sub("wikipathways\\-\\d+\\-gmt\\-([_A-Za-z]+)\\.gmt", "\\1", gmtfile)
    sub("_", " ",  orgs)
}

##' @importFrom gson read.gmt.wp
get_wp_data <- function(organism, output = "data.frame") {
    organism <- sub(" ", "_", organism)
    gmtfile <- get_wp_gmtfile()
    wpurl <- 'https://data.wikipathways.org/current/gmt/'
    url <- paste0(wpurl,
                  gmtfile[grep(organism, gmtfile)])
    f <- tempfile(fileext = ".gmt")
    dl <- mydownload(url, destfile = f)
    if (is.null(f)) {
        message("fail to download wikiPathways data...")
        return(NULL)
    }
    read.gmt.wp(f, output = output)
}
