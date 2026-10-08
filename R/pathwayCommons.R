#' ORA analysis for Pathway Commons
#'
#' This function performs over-representation analysis using  Pathway Commons
#' @title enrichPC
#' @param gene a vector of genes (either hgnc symbols or uniprot IDs)
## @param source Data source of Pathway Commons, e.g., 'reactome', 'kegg', 'pathbank', 'netpath', 'panther', etc.
## @param keyType specify the type of input 'gene' (one of 'hgnc' or 'uniprot')
#' @param ... additional parameters, see also the parameters supported by the enricher() function
#' @return A \code{enrichResult} instance
#' @importFrom utils stack
#' @export
enrichPC <- function(gene, ...) {
    # keyType <- match.arg(keyType, c("hgnc", "uniprot"))
    # source <- match.arg(source, get_pc_source())
    # pcdata <- get_pc_data(source, keyType, output = 'gson')

    pcdata <- get_pc_data(output = 'gson')
    res <- enricher(gene, gson = pcdata, ...)

    if (is.null(res)) {
        return(res)
    }

    res@ontology <- pcdata@gsname
    res@organism <- pcdata@species
    # res@keytype <-  keyType

    return(res)
}

#' GSEA analysis for  Pathway Commons
#'
#' This function performs GSEA using  Pathway Commons
#' @title gsePC
#' @param geneList a ranked gene list
## @param source Data source of Pathway Commons, e.g., 'reactome', 'kegg', 'pathbank', 'netpath', 'panther', etc.
## @param keyType specify the type of input 'gene' (one of 'hgnc' or 'uniprot')
#' @param eps boundary for calculating the p value in multilevel mode
#' @param ... additional parameters, see also the parameters supported by the GSEA() function
#' @importFrom rlang check_installed
#' @return A \code{gseaResult} instance
#' @export
gsePC <- function(geneList, eps = 1e-10, ...) {
    # keyType <- match.arg(keyType, c("hgnc", "uniprot"))
    # source <- match.arg(source, get_pc_source())

    # pcdata <- get_pc_data(source, keyType, output = 'gson')
    pcdata <- get_pc_data(output = 'gson')
    res <- GSEA(geneList, gson = pcdata, eps = eps, ...)

    if (is.null(res)) {
        return(res)
    }

    res@ontology <- pcdata@gsname
    res@organism <- pcdata@species
    # res@keytype <-  keyType

    return(res)
}

prepare_pc_data <- function() {
    pc2gene <- get_pc_data(output = 'data.frame')
    ##TERM2GENE
    pcid2gene <- pc2gene[, c("id", "gene")]
    ##TERM2NAME
    pcid2name <- unique(pc2gene[, c("id", "name")])

    list(PCID2GENE = pcid2gene, PCID2NAME = pcid2name)
}

get_pc_gmtfile <- function() {
    # pcurl <- 'https://www.pathwaycommons.org/archives/PC2/v12/'
    pc2 <- 'https://download.baderlab.org/PathwayCommons/PC2/'
    con <- readLines(pc2)
    pattern <- ".*>v(\\d+)/</a>.*"
    latest_version <- con[grep(pattern, con)] |>
        sub(pattern, "\\1", x = _) |>
        as.numeric() |>
        max()

    pcurl <- sprintf('%sv%s/', pc2, latest_version)

    x <- readLines(pcurl)
    y <- x[grep('\\.gmt.gz', x)]
    # sub(".*(PathwayCommons.*\\.gmt.gz).*", "\\1",  y)
    file <- sub(".*>(.*\\.gmt\\.gz)</a>.*", "\\1", y)
    sprintf("%s%s", pcurl, file)
}

### list supported data sources of Pathway Commons
# get_pc_source <- function() {
#     gmtfile <- get_pc_gmtfile()
#     source <- unique(sub("PathwayCommons\\d+\\.([_A-Za-z]+)\\.([_A-Za-z]+)\\.gmt.gz", "\\1", gmtfile))
#
#     return(source)
# }

read.gmt.pc_internal <- function(gmtfile) {
    check_installed(
        'readr',
        'for `read.gmt.pc_internal()`, which is an internal function.'
    )

    x <- yread(gmtfile, readr::read_lines)
    y <- strsplit(x, "\t")

    url <- vapply(y, `[`, 1, FUN.VALUE = character(1))
    description <- vapply(y, `[`, 2, FUN.VALUE = character(1))
    ## a blank line holds no gene set at all, and must not carry a negative count
    ngene <- pmax(vapply(y, length, FUN.VALUE = integer(1)) - 2L, 0L)

    data.frame(
        id = rep(sub(".*/", "", url), ngene),
        description = rep(description, ngene),
        gene = unlist(lapply(y, `[`, -c(1:2)), use.names = FALSE),
        stringsAsFactors = FALSE
    )
}

#' The labelled fields of a Pathway Commons gene set description, which reads
#' `name: ...; datasource: ...; organism: ...; idtype: ...`
#'
#' @noRd
pc_field <- function(description, key) {
    label <- regexpr(paste0("(^|; )", key, ": "), description)
    value <- sub(";.*$", "", substring(description, label + attr(label, "match.length")))
    ifelse(label > 0, value, NA_character_)
}

#' @noRd
pc_name <- function(description) {
    ## only the name can hold "; " itself, so it ends at a known label, not at the first semicolon
    name <- sub("^name: ", "", description)
    sub("; (datasource|organism|idtype): .*", "", name)
}

#' Parse gmt file from Pathway Common
#'
#' This function parse gmt file downloaded from Pathway common
#' @title read.gmt.pc
#' @param gmtfile A gmt file
#' @param output one of 'data.frame' or 'GSON'
#' @return A data.frame or A GSON object depends on the value of 'output'
#' @export
read.gmt.pc <- function(gmtfile, output = "data.frame") {
    output <- match.arg(output, c("data.frame", "gson", "GSON"))

    pcdata <- read.gmt.pc_internal(gmtfile)
    x <- data.frame(
        id = pcdata$id,
        name = pc_name(pcdata$description),
        datasource = pc_field(pcdata$description, "datasource"),
        organism = pc_field(pcdata$description, "organism"),
        idtype = pc_field(pcdata$description, "idtype"),
        gene = pcdata$gene,
        stringsAsFactors = FALSE
    )

    if (output == "data.frame") {
        return(x)
    }

    gsid2gene <- data.frame(gsid = x$id, gene = x$gene)
    gsid2name <- unique(data.frame(gsid = x$id, name = x$name))
    taxid <- x$organism[!is.na(x$organism)]
    organism <- if (length(taxid) > 0) taxID2name(taxid[1]) else NA_character_
    gson(
        gsid2gene = gsid2gene,
        gsid2name = gsid2name,
        gsname = "Pathway Commons",
        species = organism
    )
}


get_pc_data <- function(output = "data.frame") {
    url <- get_pc_gmtfile()
    # gmtfile <- gmtfile[grepl(source, gmtfile) & grepl(keyType, gmtfile)]

    read.gmt.pc(url, output = output)
}
