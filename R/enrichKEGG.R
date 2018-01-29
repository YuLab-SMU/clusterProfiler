##' KEGG Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment KEGG
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff
##' @param use_internal_data logical, use KEGG.db or latest online KEGG data
##' @return A \code{enrichResult} instance.
##' @export
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @importClassesFrom methods data.frame
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @examples
##'
##' 	data(geneList, package='DOSE')
##'   de <- names(geneList)[1:100]
##' 	yy <- enrichKEGG(de, pvalueCutoff=0.01)
##' 	head(yy)
##'
##'
enrichKEGG <- function(gene,
                       organism          = "hsa",
                       keyType           = "kegg",
                       pvalueCutoff      = 0.05,
                       pAdjustMethod     = "BH",
                       universe,
                       minGSSize         = 10,
                       maxGSSize         = 500,
                       qvalueCutoff      = 0.2,
                       use_internal_data = FALSE) {

    species <- organismMapper(organism)
    if (use_internal_data) {
        KEGG_DATA <- get_data_from_KEGG_db(species)
    } else {
        KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
    }
    res <- enricher_internal(gene,
                             pvalueCutoff  = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             maxGSSize     = maxGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             USER_DATA = KEGG_DATA)
    if (is.null(res))
        return(res)

    res@ontology <- "KEGG"
    res@organism <- species
    res@keytype <- keyType

    return(res)
}

get_KEGG_Env <- function() {
    if (! exists(".KEGG_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos)
        assign(".KEGG_clusterProfiler_Env", new.env(), envir=envir)
    }
    get(".KEGG_clusterProfiler_Env", envir = .GlobalEnv)
}

##' download the latest version of KEGG pathway/module
##'
##'
##' @title download_KEGG
##' @param species species
##' @param keggType one of 'KEGG' or 'MKEGG'
##' @param keyType supported keyType, see bitr_kegg
##' @return list
##' @author Guangchuang Yu
##' @importFrom magrittr %<>%
##' @export
download_KEGG <- function(species, keggType="KEGG", keyType="kegg") {
    KEGG_Env <- get_KEGG_Env()

    use_cached <- FALSE

    if (exists("organism", envir = KEGG_Env, inherits = FALSE) &&
        exists("_type_", envir = KEGG_Env, inherits = FALSE) ) {

        org <- get("organism", envir=KEGG_Env)
        type <- get("_type_", envir=KEGG_Env)

        if (org == species && type == keggType &&
            exists("KEGGPATHID2NAME", envir=KEGG_Env, inherits = FALSE) &&
            exists("KEGGPATHID2EXTID", envir=KEGG_Env, inherits = FALSE)) {

            use_cached <- TRUE
        }
    }

    if (use_cached) {
        KEGGPATHID2EXTID <- get("KEGGPATHID2EXTID", envir=KEGG_Env)
        KEGGPATHID2NAME <- get("KEGGPATHID2NAME", envir=KEGG_Env)
    } else {
        if (keggType == "KEGG") {
            kres <- download.KEGG.Path(species)
        } else {
            kres <- download.KEGG.Module(species)
        }

        KEGGPATHID2EXTID <- kres$KEGGPATHID2EXTID
        KEGGPATHID2NAME <- kres$KEGGPATHID2NAME

        assign("organism", species, envir=KEGG_Env)
        assign("_type_", keggType, envir=KEGG_Env)
        assign("KEGGPATHID2NAME", KEGGPATHID2NAME, envir=KEGG_Env)
        assign("KEGGPATHID2EXTID", KEGGPATHID2EXTID, envir=KEGG_Env)
    }

    if (keyType != "kegg") {
        need_idconv <- FALSE
        idconv <- NULL
        if (use_cached &&
            exists("key", envir=KEGG_Env, inherits = FALSE) &&
            exists("idconv", envir=KEGG_Env, inherits = FALSE)) {

            key <- get("key", envir=KEGG_Env)
            if (key == keyType) {
                idconv <- get("idconv", envir=KEGG_Env)
            } else {
                need_idconv <- TRUE
            }
        } else {
            neec_idconv <- TRUE
        }

        if (need_idconv || is.null(idconv)) {
            idconv <- KEGG_convert("kegg", keyType, species)
            assign("key", keyType, envir=KEGG_Env)
            assign("idconv", idconv, envir=KEGG_Env)
        }
        colnames(KEGGPATHID2EXTID) <- c("from", "kegg")
        KEGGPATHID2EXTID <- merge(KEGGPATHID2EXTID, idconv, by.x='kegg', by.y='from')
        KEGGPATHID2EXTID <- unique(KEGGPATHID2EXTID[, -1])
    }

    return(list(KEGGPATHID2EXTID = KEGGPATHID2EXTID,
                KEGGPATHID2NAME  = KEGGPATHID2NAME))
}

prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    kegg <- download_KEGG(species, KEGG_Type, keyType)
    build_Anno(kegg$KEGGPATHID2EXTID,
               kegg$KEGGPATHID2NAME)
}

download.KEGG.Path <- function(species) {
    keggpathid2extid.df <- kegg_link(species, "pathway")
    if (is.null(keggpathid2extid.df))
        stop("'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...")
    keggpathid2extid.df[,1] %<>% gsub("[^:]+:", "", .)
    keggpathid2extid.df[,2] %<>% gsub("[^:]+:", "", .)

    keggpathid2name.df <- kegg_list("pathway")
    keggpathid2name.df[,1] %<>% gsub("path:map", species, .)

    ## if 'species="ko"', ko and map path are duplicated, only keep ko path.
    ##
    ## http://www.kegg.jp/dbget-bin/www_bget?ko+ko00010
    ## http://www.kegg.jp/dbget-bin/www_bget?ko+map0001
    ##
    keggpathid2extid.df <- keggpathid2extid.df[keggpathid2extid.df[,1] %in% keggpathid2name.df[,1],]

    return(list(KEGGPATHID2EXTID=keggpathid2extid.df,
                KEGGPATHID2NAME=keggpathid2name.df))
}

download.KEGG.Module <- function(species) {
    keggmodule2extid.df <- kegg_link(species, "module")
    if (is.null(keggmodule2extid.df)) {
        stop("'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...")
    }

    keggmodule2extid.df[,1] %<>% gsub("[^:]+:", "", .) %>% gsub(species, "", .) %>% gsub("^_", "", .)
    keggmodule2extid.df[,2] %<>% gsub("[^:]+:", "", .)

    keggmodule2name.df <- kegg_list("module")
    keggmodule2name.df[,1] %<>% gsub("md:", "", .)
    return(list(KEGGPATHID2EXTID=keggmodule2extid.df,
                KEGGPATHID2NAME =keggmodule2name.df))
}


##' viewKEGG function is for visualize KEGG pathways
##' works with enrichResult object to visualize enriched KEGG pathway
##'
##'
##' @param obj enrichResult object
##' @param pathwayID pathway ID or index
##' @param foldChange fold change values
##' @param color.low color of low foldChange genes
##' @param color.high color of high foldChange genes
##' @param kegg.native logical
##' @param out.suffix suffix of output file
## @importFrom pathview pathview
## @importFrom pathview kegg.species.code
##' @importFrom utils citation
##' @references Luo et al. (2013) Pathview: an R/Bioconductor package for
##'pathway-based data integration and visualization. \emph{Bioinformatics} (Oxford,
##'England), 29:14 1830--1831, 2013. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/content/abstract/29/14/1830.abstract}
##'PMID: 23740750
viewKEGG <- function(obj, pathwayID, foldChange,
                       color.low="green",
                       color.high="red",
                       kegg.native=TRUE,
                       out.suffix="clusterProfiler") {

    if (class(obj) != "enrichResult")
        stop("only enrichResult object supported.")
    if (obj@ontology != "KEGG")
        stop("only KEGG supported.")

    print("viewKEGG is a wrapper function of pathview")
    citation("pathview")

    pkg <- "pathview"
    suppressMessages(require(pkg, character.only=TRUE))
    if (is.numeric(pathwayID)) {
        pathwayID <- as.data.frame(obj)[pathwayID, 1]
    }
    if (length(pathwayID) == 1 & pathwayID == "all") {
        pathwayID <- as.data.frame(obj)[, 1]
    }
    m.fc <- max(abs(foldChange))
    bins <- ceiling(m.fc) * 2
    if (bins < 10)
        bins <- 10
    pathview <- eval(parse(text=pkg))
    res <- lapply(pathwayID, function(pid) {
        pathview(gene.data=foldChange,
                 pathway.id = pid,
                 species = "hsa",
                 limit = list(gene=m.fc, cpd=1),
                 bins = list(gene=bins, cpd=10),
                 low = list(gene=color.low, cpd="blue"),
                 high = list(gene=color.high, cpd="yellow"),
                 kegg.native=kegg.native,
                 out.suffix=out.suffix,
                 new.signature=FALSE)
    })
    return (res)
}

##' @importFrom AnnotationDbi as.list
##' @importFrom utils stack
get_data_from_KEGG_db <- function(species) {
    PATHID2EXTID <- as.list(get_KEGG_db("KEGGPATHID2EXTID"))
    if (!any(grepl(species, names(PATHID2EXTID)))) {
        stop("input species is not supported by KEGG.db...")
    }
    idx <- grep(species, names(PATHID2EXTID))
    PATHID2EXTID <- PATHID2EXTID[idx]
    PATHID2EXTID.df <- stack(PATHID2EXTID)
    PATHID2EXTID.df <- PATHID2EXTID.df[, c(2,1)]
    PATHID2NAME <- as.list(get_KEGG_db("KEGGPATHID2NAME"))
    PATHID2NAME.df <- data.frame(path=paste0(species, names(PATHID2NAME)),
                                 name=unlist(PATHID2NAME))
    build_Anno(PATHID2EXTID.df, PATHID2NAME.df)
}

get_KEGG_db <- function(kw) {
    annoDb <- "KEGG.db"
    suppressMessages(requireNamespace(annoDb))
    eval(parse(text=paste0(annoDb, "::", kw)))
}

organismMapper <- function(organism) {
    ## it only map those previous supported organism

    if (organism == "anopheles") {
        species <- "aga"
    } else if (organism == "arabidopsis") {
        species <- "ath"
    } else if (organism == "bovine") {
        species <- "bta"
    } else if (organism == "canine") {
        species <- "cfa"
    } else if (organism == "chicken") {
        species <- "gga"
    } else if (organism == "chipm") {
        species <- "ptr"
    } else if (organism == "ecolik12") {
        species <- "eco"
    } else if (organism == "ecsakai") {
        species <- "ecs"
    } else if (organism == "fly") {
        species <- "dme"
    } else if (organism == "human") {
        species <- "hsa"
    } else if (organism == "malaria") {
        species <- "pfa"
    } else if (organism == "mouse") {
        species <- "mmu"
    } else if (organism == "pig") {
        species <- "ssc"
    } else if (organism == "rat") {
        species <- "rno"
    } else if (organism == "rhesus") {
        species <- "mcc"
    } else if (organism == "worm" || organism == "celegans") {
        species <- "cel"
    } else if (organism == "xenopus") {
        species <- "xla"
    } else if (organism == "yeast") {
        species <- "sce"
    } else if (organism == "zebrafish") {
        species <- "dre"
    } else {
        species <- organism
    }
    return(species)
}




