##' download the latest version of KEGG pathway and stored in a 'GSON' object
##'
##'
##' @title gson_KEGG
##' @param species species
##' @param KEGG_Type one of "KEGG" and "MKEGG"
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'.
##' @return a 'GSON' object
##' @author Guangchuang Yu
##' @importFrom gson gson
##' @export
gson_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    x <- download_KEGG(species, KEGG_Type, keyType)
    gsid2gene <- setNames(x[[1]], c("gsid", "gene"))
    gsid2name <- setNames(x[[2]], c("gsid", "name"))
    y <- readLines("https://rest.kegg.jp/info/hsa")
    version <- sub("\\w+\\s+", "", y[grep('Release', y)])
    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name,
        species = species,
        gsname = "KEGG",
        version = version,
        accessed_date = as.character(Sys.Date())
    )
}


gson_GO <- function(OrgDb, keytype = 'ENTREZID', ont = "BP") {

    if (is(OrgDb, "character")) {
        require(OrgDb, character.only = TRUE)
        OrgDb <- eval(parse(text = OrgDb))
    }

    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    if (ont != "ALL") {
        goterms <- goterms[goterms == ont]
    }
    go2gene <- suppressMessages(
        AnnotationDbi::mapIds(OrgDb, keys=names(goterms), column=keytype,
                            keytype="GOALL", multiVals='list')
    )
    goAnno <- stack(go2gene)   
    gsid2gene <- goAnno[, c(2,1)]
    colnames(gsid2gene) <- c("gsid", "gene")
    gsid2gene <- unique(gsid2gene[!is.na(gsid2gene[,2]), ]) 

    termname <- AnnotationDbi::Term(GO.db::GOTERM)
    gsid2name <- data.frame(gsid = names(termname),
                            name = termname)
    species <- AnnotationDbi::species(OrgDb)
    m <- AnnotationDbi::metadata(OrgDb)
    version <- m$value[m$name == "GOSOURCEDATE"]
    # gsname <- m$value[m$name == 'GOSOURCENAME']
    gsname <- paste(m$value[m$name == 'GOSOURCENAME'], ont, sep = ";")
    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name,
        species = species,
        gsname = gsname,
        version = version,
        accessed_date = as.character(Sys.Date()),
        keytype = keytype
    )
}
##' Download the latest version of WikiPathways data and stored in a 'GSON' object
##'
##'
##' @title gson_WP
##' @param organism supported organism, which can be accessed via the get_wp_organisms() function.
##' @export
gson_WP <- function(organism) {
    get_wp_data(organism, output = "GSON")
}
