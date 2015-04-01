.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
}

build_Anno <- function(path2gene, path2name) {
    if (!exists("Anno_clusterProfiler_Env", envir = .GlobalEnv)) {
        assign("Anno_clusterProfiler_Env", new.env(), .GlobalEnv)
    }
    Anno_clusterProfiler_Env <- get("Anno_clusterProfiler_Env", envir= .GlobalEnv)

    PATHID2EXTID <- split(as.character(path2gene[,2]), as.character(path2gene[,1]))
    EXTID2PATHID <- split(as.character(path2gene[,1]), as.character(path2gene[,2]))
    
    assign("PATHID2EXTID", PATHID2EXTID, envir = Anno_clusterProfiler_Env)
    assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)

    if (! missing(path2name) || is.null(path2name) || is.na(path2name)) {
        assign("PATHID2NAME", path2name, envir = Anno_clusterProfiler_Env)
    } else {
        assign("PATHID2NAME", NULL, envir = Anno_clusterProfiler_Env)
    }
    return(Anno_clusterProfiler_Env)
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

KEGG_db_supported <- function() {
    res <- c(
        "aga",
        "ath",
        "bta",
        "cfa",
        "gga",
        "ptr",
        "eco",
        "ecs",
        "dme",
        "hsa",
        "pfa",
        "mmu",
        "ssc",
        "rno",
        "mcc",
        "cel",
        "xla",
        "sce",
        "dre")

    return(res)
}



get_KEGG_Anno <- function(organism, key) {
    buildKEGGmap_online(organism)
    KEGG_clusterProfiler_Env <- get("KEGG_clusterProfiler_Env", envir = .GlobalEnv)
    get(key, KEGG_clusterProfiler_Env)
}

buildKEGGmap_online <- function(organism) {
    organism <- organismMapper(organism)

    if (! exists("KEGG_clusterProfiler_Env", envir = .GlobalEnv)) {
        assign("KEGG_clusterProfiler_Env", new.env(), .GlobalEnv)
    }
    KEGG_clusterProfiler_Env <- get("KEGG_clusterProfiler_Env", envir = .GlobalEnv)
    
    if (exists("organism", envir = KEGG_clusterProfiler_Env, inherits = FALSE)) {
        org <- get("organism", envir=KEGG_clusterProfiler_Env)
        if (org == organism &&
            exists("KEGGPATHID2NAME", envir=KEGG_clusterProfiler_Env, inherits = FALSE) &&
            exists("KEGGPATHID2EXTID", envir=KEGG_clusterProfiler_Env, inherits = FALSE) &&
            exists("EXTID2KEGGPATHID", envir=KEGG_clusterProfiler_Env, inherits = FALSE)
            ) {
            ## do nothing
        } else {
            download_buildKEGGmap(organism)
        }
    } else {
        download_buildKEGGmap(organism)
    }
    
}


download_buildKEGGmap <- function(organism) {
    organism <- organismMapper(organism)
    kegg <- download.KEGG(organism)
    buildKEGGmap(kegg[[1]], kegg[[2]], organism)
}


##' download the latest version of KEGG pathway
##'
##' 
##' @title download.KEGG
##' @param species species
##' @return list
##' @author Guangchuang Yu
##' @importFrom KEGGREST keggLink
##' @importFrom KEGGREST keggList
##' @importFrom magrittr %<>%
download.KEGG <- function(species) {
    keggpathid2extid <- keggLink(species,"pathway")
    keggpathid2extid %<>% gsub("[^:]+:", "", .)
    names(keggpathid2extid) %<>% gsub("[^:]+:", "", .)

    keggpath2extid.df <- data.frame(pathID=names(keggpathid2extid), extID=keggpathid2extid)
    
    keggpathid2name<-keggList("pathway")
    names(keggpathid2name) %<>% gsub("path:map", "", .)

    res <- list(keggpath2extid = keggpath2extid.df,
                keggpathid2name = keggpathid2name)
    return(res)
}


##' build KEGG annotation files
##'
##' 
##' @title buildKEGGmap
##' @param keggmap pathway to external ID
##' @param id2name pathway id to pathway name
##' @param organism organism
##' @return NULL
##' @importFrom magrittr %>%
##' @author Guangchuang Yu
buildKEGGmap <- function(keggmap, id2name=NULL, organism) {
    if (! exists("KEGG_clusterProfiler_Env", envir = .GlobalEnv)) {
        assign("KEGG_clusterProfiler_Env", new.env(), .GlobalEnv)
    }
    KEGG_clusterProfiler_Env <- get("KEGG_clusterProfiler_Env", envir = .GlobalEnv)

    if (is.null(id2name)) {
        id2name <- as.list(KEGGPATHID2NAME) %>% unlist
        pathid  <- names(id2name)
        id      <- keggmap[,1] %>% as.character %>% gsub("^[a-z]+", "", .)
        keggmap <- keggmap[id %in% pathid, ]
        
    }
    
    KEGGPATHID2NAME  <- id2name
    KEGGPATHID2EXTID <- split(as.character(keggmap[,2]), as.character(keggmap[,1]))
    EXTID2KEGGPATHID <- split(as.character(keggmap[,1]), as.character(keggmap[,2]))

    assign("organism", organism, envir = KEGG_clusterProfiler_Env)
    assign("KEGGPATHID2NAME", KEGGPATHID2NAME,  envir  = KEGG_clusterProfiler_Env)
    assign("KEGGPATHID2EXTID", KEGGPATHID2EXTID, envir = KEGG_clusterProfiler_Env)
    assign("EXTID2KEGGPATHID", EXTID2KEGGPATHID, envir = KEGG_clusterProfiler_Env)
}


excludeGOlevel <- function(x, ont, level) {
    lv <- unlist(lapply(level, getGOLevel, ont=ont))
    x <- excludeGOterm(x, lv)
    return(x)
}

excludeGOterm <- function(x, term) {
    if ( is(x, "enrichResult") ) {
        x@result <- x@result[! x@result[, "ID"] %in% term, ]
    } else if ( is(x, "compareClusterResult") ) {
        x@compareClusterResult <- x@compareClusterResult[! x@compareClusterResult[, "ID"] %in% term, ]
    } else {
        stop("x should be one of enrichResult of compareClusterResult...")
    }
    return(x)
}

keepGOlevel <- function(x, ont, level) {
    lv <- unlist(lapply(level, getGOLevel, ont=ont))
    x <- keepGOterm(x, lv)
    return(x)
}

keepGOterm <- function(x, term) {
    if ( is(x, "enrichResult") ) {
        x@result <- x@result[x@result[, "ID"] %in% term, ]
    } else if ( is(x, "compareClusterResult") ) {
        x@compareClusterResult <- x@compareClusterResult[x@compareClusterResult[, "ID"] %in% term, ]
    } else {
        stop("x should be one of enrichResult of compareClusterResult...")
    }
    return(x)
}

##' @importFrom GOSemSim getDb
getGO2ALLEG_MappedDb <- function(organism) {
    annoDb <- getDb(organism)
    require(annoDb, character.only = TRUE)

    mappedDb <- switch(organism,
                       anopheles    = "org.Ag.egGO2ALLEGS",
                       arabidopsis  = "org.At.tairGO2ALLTAIRS",
                       bovine       = "org.Bt.egGO2ALLEGS",
                       canine       = "org.Cf.egGO2ALLEGS",
                       chicken      = "org.Gg.egGO2ALLEGS",
                       chimp        = "org.Pt.egGO2ALLEGS",
                       coelicolor   = "org.Sco.egGO2ALLEGS",
                       ecolik12     = "org.EcK12.egGO2ALLEGS",
                       ecsakai      = "org.EcSakai.egGO2ALLEGS",
                       fly          = "org.Dm.egGO2ALLEGS",
                       gondii       = "org.Tgondii.egGO2ALLEGS",
                       human        = "org.Hs.egGO2ALLEGS",
                       malaria      = "org.Pf.plasmoGO2ALLORFS",
                       mouse        = "org.Mm.egGO2ALLEGS",
                       pig          = "org.Ss.egGO2ALLEGS",
                       rat          = "org.Rn.egGO2ALLEGS",
                       rhesus       = "org.Mmu.egGO2ALLEGS",
                       worm         = "org.Ce.egGO2ALLEGS",
                       xenopus      = "org.Xl.egGO2ALLEGS",
                       yeast        = "org.Sc.sgdGO2ALLORFS",
                       zebrafish    = "org.Dr.egGO2ALLEGS"
                       )
    mappedDb <- eval(parse(text=mappedDb))
    return(mappedDb)
}

##' @importFrom GOSemSim loadGOMap
getEG2GO_MappedDb <- function(organism) {
    loadGOMap(organism)
    mappedDb <- get("gomap", envir=eval(parse(text="GOSemSimEnv")))
    return(mappedDb)
}


##' query GOIDs at a specific level.
##'
##'
##' @title get GOIDs at a specific level
##' @param ont Ontology
##' @param level GO level
##' @return a vector of GOIDs
##' @importFrom GO.db GOBPCHILDREN
##' @importFrom GO.db GOCCCHILDREN
##' @importFrom GO.db GOMFCHILDREN
##' @importMethodsFrom AnnotationDbi mget
##' @author Guangchuang Yu \url{http://ygc.name}
getGOLevel <- function(ont, level) {
    switch(ont,
           MF = {
               topNode <- "GO:0003674"
               Children <- GOMFCHILDREN
           },
           BP = {
               topNode <- "GO:0008150"
               Children <- GOBPCHILDREN
           },
           CC = {
               topNode <- "GO:0005575"
               Children <- GOCCCHILDREN
           }
           )

    Node <- topNode
    for (i in seq_len(level-1)) {
        Node <- mget(Node, Children, ifnotfound=NA)
        Node <- unique(unlist(Node))
        Node <- as.vector(Node)
        Node <- Node[!is.na(Node)]
    }
    return(Node)
}


##' Internal plot function for plotting compareClusterResult
##'
##'
##' @title plotting-clusterProfile
##' @param clProf.reshape.df data frame of compareCluster result
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @param colorBy one of pvalue or p.adjust
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 scale_colour_gradient
##' @importFrom DOSE theme_dose
##' @author Guangchuang Yu \url{http://ygc.name}
plotting.clusterProfile <- function(clProf.reshape.df,
                                    type = "dot",
                                    by = "geneRatio",
                                    colorBy = "p.adjust",
                                    title="",
                                    font.size=12) {
    Description <- Percentage <- Count <- Cluster <- GeneRatio <- p.adjust <- pvalue <- NULL # to satisfy codetools
    if (type == "bar") {
        if (by == "percentage") {
            p <- ggplot(clProf.reshape.df,
                        aes(x=Description, y = Percentage, fill=Cluster))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df,
                        aes(x=Description, y = Count, fill=Cluster))
        } else {

        }
        p <- p +
            geom_bar() +
                coord_flip()
    }
    if (type == "dot") {
        if (by == "rowPercentage") {
            p <- ggplot(clProf.reshape.df,
                        aes(x = Cluster, y = Description, size = Percentage))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df,
                        aes(x = Cluster, y = Description, size = Count))
        } else if (by == "geneRatio") {
            p <- ggplot(clProf.reshape.df,
                        aes(x = Cluster, y = Description, size = GeneRatio))
        } else {
            ## nothing here
        }
        if (any(colnames(clProf.reshape.df) == colorBy)) {
            p <- p +
                geom_point() +
                    aes_string(color=colorBy) +
                        scale_colour_gradient(low="red", high="blue")
        } else {
            p <- p + geom_point(colour="steelblue")
        }
    }
    p <- p + xlab("") + ylab("") + ggtitle(title) +
        theme_dose(font.size)
    ## theme(axis.text.x = element_text(colour="black", size=font.size, vjust = 1)) +
    ##     theme(axis.text.y = element_text(colour="black",
    ##           size=font.size, hjust = 1)) +
    ##               ggtitle(title)+theme_bw()
    ## p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
    ##                    hjust=hjust.axis.x,
    ##                    vjust=vjust.axis.x))
    return(p)
}

##' building GO mapping files
##'
##' provided by a data.frame of gene and GO directly annotation file
##' this function will building gene to GO and GO to gene mapping,
##' with directly and undirectly annotation.
##' @title buildGOmap
##' @param gomap data.frame with two columns names "entrezgene", and "go_accession"
##' @param compress logical, indicate file save in compress or not.
##' @return files save in the the working directory
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
##' @importFrom plyr dlply
##' @export
##' @author Yu Guangchuang
buildGOmap <- function(gomap, compress=TRUE) {
    if( any( colnames(gomap) %in% "go_id" ) ) {
        colnames(gomap)[colnames(gomap) %in% "go_id"] <- "go_accession"
    }

    ## remove empty GO annotation
    gomap <- gomap[gomap$go_accession != "", ]


    GO2EG <- dlply(gomap,"go_accession",.fun=function(i) as.character(i$entrezgene))
    EG2GO <- dlply(gomap,"entrezgene",.fun=function(i) as.character(i$go_accession))

    if (compress) {
      save(GO2EG, file="GO2EG.rda", compress="xz")
      save(EG2GO, file="EG2GO.rda", compress="xz")
    } else {
      save(GO2EG, file="GO2EG.rda")
      save(EG2GO, file="EG2GO.rda")

    }


    EG2ALLGO <- lapply(EG2GO,
                        function(i) {
                            mfans <- unlist(mget(i, GOMFANCESTOR, ifnotfound=NA))
                            bpans <- unlist(mget(i, GOBPANCESTOR, ifnotfound=NA))
                            ccans <- unlist(mget(i, GOCCANCESTOR, ifnotfound=NA))
                            ans <- c(mfans, bpans, ccans)
                            ans <- ans[ !is.na(ans) ]
                            ans <- c(i, ans)
                            ans <- unique(ans)
                            ans <- ans[ans != "all"]
                            return(ans)
                        })
    if (compress) {
      save(EG2ALLGO, file="EG2ALLGO.rda", compress="xz")
    } else {
      save(EG2ALLGO, file="EG2ALLGO.rda")
    }

    len <- lapply(EG2ALLGO,length)
    EG2ALLGO.df <- data.frame(EG=rep(names(EG2ALLGO), times=len),
                              GO=unlist(EG2ALLGO))
    GO <- NULL ## satisfy code tools
    GO2ALLEG <- dlply(EG2ALLGO.df, .(GO), function(i) as.character(i$EG))
    GO2ALLEG <- lapply(GO2ALLEG, unique)
    if (compress) {
      save(GO2ALLEG, file="GO2ALLEG.rda", compress="xz")
    } else {
      save(GO2ALLEG, file="GO2ALLEG.rda")
    }
    print("GO Annotation Mapping files save in the working directory.")
}

GI2EG <- function(GI, organism="D39") {
    gi <- as.character(GI)
    ## remove blank
    gi <- sub("^\\s+", "", gi, perl=T)
    gi <- sub("\\s+$", "", gi, perl=T)
    ## remove GI: or gi:
    gi <- sapply(gi, function(i) unlist(strsplit(i, split="\\|"))[2])
    ## load corresponding gene table and protein table
    geneTable <- proteinTable <- NULL
    if (organism == "M5005") {
        f <- system.file("extdata", "M5005/geneTable.rda", package="clusterProfiler")
        load(f)
        gi.eg <- geneTable[geneTable$GI %in% gi, c("GI", "GeneID")]
    } else if (organism=="D39") {
        gt <- system.file("extdata", "D39/geneTable.rda", package="clusterProfiler")
        load(gt)
        pt <- system.file("extdata", "D39/proteinTable.rda", package="clusterProfiler")
        load(pt)
        idx <- match(gi, proteinTable$PID)
        locus <- proteinTable[idx, "Synonym"]
        locus <- as.character(locus)
        idx <- match(locus, geneTable$Locus)
        gene <- geneTable[idx, "GeneID"]
        gene <- as.character(gene)
        gi.eg <- data.frame(GI=gi, GeneID=gene)
    } else {
        stop("not supported yet...")
    }
    
    return(gi.eg)
}


removeEmptyEntry.list <- function(x) {
    notNA.idx <- unlist(lapply(x, function(i) !is.null(i) && !all(is.na(i))))
    x[notNA.idx]
}
