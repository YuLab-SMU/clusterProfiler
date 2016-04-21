.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
}

build_Anno <- function(path2gene, path2name) {
    if (!exists(".Anno_clusterProfiler_Env", envir = .GlobalEnv)) {
        assign(".Anno_clusterProfiler_Env", new.env(), .GlobalEnv)
    }
    Anno_clusterProfiler_Env <- get(".Anno_clusterProfiler_Env", envir= .GlobalEnv)

    PATHID2EXTID <- split(as.character(path2gene[,2]), as.character(path2gene[,1]))
    EXTID2PATHID <- split(as.character(path2gene[,1]), as.character(path2gene[,2]))
    
    assign("PATHID2EXTID", PATHID2EXTID, envir = Anno_clusterProfiler_Env)
    assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)

    if ( missing(path2name) || is.null(path2name) || is.na(path2name)) {
        assign("PATHID2NAME", NULL, envir = Anno_clusterProfiler_Env)
    } else {
	path2name <- unique(path2name)
	PATH2NAME <- as.character(path2name[,2])
	names(PATH2NAME) <- as.character(path2name[,1]) 
        assign("PATHID2NAME", PATH2NAME, envir = Anno_clusterProfiler_Env)
    }
    return(Anno_clusterProfiler_Env)
}



## buildKEGGmap <- function(keggmap, id2name=NULL, organism) {
##     if (! exists("KEGG_clusterProfiler_Env", envir = .GlobalEnv)) {
##         assign("KEGG_clusterProfiler_Env", new.env(), .GlobalEnv)
##     }
##     KEGG_clusterProfiler_Env <- get("KEGG_clusterProfiler_Env", envir = .GlobalEnv)

##     if (is.null(id2name)) {
##         id2name <- as.list(KEGGPATHID2NAME) %>% unlist
##         pathid  <- names(id2name)
##         id      <- keggmap[,1] %>% as.character %>% gsub("^[a-z]+", "", .)
##         keggmap <- keggmap[id %in% pathid, ]
        
##     }
    
##     KEGGPATHID2NAME  <- id2name
##     KEGGPATHID2EXTID <- split(as.character(keggmap[,2]), as.character(keggmap[,1]))
##     EXTID2KEGGPATHID <- split(as.character(keggmap[,1]), as.character(keggmap[,2]))

##     assign("organism", organism, envir = KEGG_clusterProfiler_Env)
##     assign("KEGGPATHID2NAME", KEGGPATHID2NAME,  envir  = KEGG_clusterProfiler_Env)
##     assign("KEGGPATHID2EXTID", KEGGPATHID2EXTID, envir = KEGG_clusterProfiler_Env)
##     assign("EXTID2KEGGPATHID", EXTID2KEGGPATHID, envir = KEGG_clusterProfiler_Env)
## }




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
##' @param x x variable
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @param colorBy one of pvalue or p.adjust
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_
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
                                    x = ~Cluster,
                                    type = "dot",
                                    colorBy = "p.adjust",
                                    by = "geneRatio",
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
                        aes_(x = x, y = ~Description, size = ~Percentage))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df,
                        aes_(x = x, y = ~Description, size = ~Count))
        } else if (by == "geneRatio") {
            p <- ggplot(clProf.reshape.df,
                        aes_(x = x, y = ~Description, size = ~GeneRatio))
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

get_go_ontology <- function(x) {
    if (is(x, "compareClusterResult")) {
        if (x@fun != "enrichGO") {
            stop("simplify only work for GO...")
        }
        ont <- x@.call$ont
        if (is.null(ont)) {
            ## should be "MF", default value of enrichGO
            ## it's safe to determine from the output
            ont <- x@compareClusterResult$ID[1] %>% GOTERM[[.]] %>% Ontology
        }
    } else if (is(x, "enrichResult")) {
        if (!x@ontology %in% c("BP", "MF", "CC"))
            stop("ontology should be one of 'MF', 'BP', 'CC'...")

        ont <- x@ontology
    } else {
        stop("x should be an instance of 'enrichResult' or 'compareClusterResult'...")
    }

    return(ont)
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


##' @importFrom S4Vectors metadata
get_organism <- function(OrgDb) {
    OrgDb <- load_OrgDb(OrgDb)
    md <- metadata(OrgDb)
    md[md[,1] == "ORGANISM", 2]
}

add_GO_Ontology <- function(obj, GO_DATA) {
    obj@setType <- "GOALL"
    df <- obj@result
    GO2ONT <- get("GO2ONT", envir=GO_DATA)
    df <- cbind(ONTOLOGY=GO2ONT[df$ID], df)
    obj@result <- df
    return(obj)
}
