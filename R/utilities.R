
##' Assign some initial instances.
##'
##'
##' @title initial function for clusterProfiler
##' @return environment instance
##' @author Guangchuang Yu
.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
}

##' provide a vector of GOIDs, this function will convert them to corresponding GO Terms
##'
##'
##' @title Mapping GOIDs to GO Terms
##' @param GOID GOID
##' @return GO Terms
##' @author Guangchuang Yu
GO2Term <- function(GOID) {
    go <- mget(GOID, GOTERM, ifnotfound=NA)
    term <- sapply(go, Term)
    return(term)
}


##' provide a vector of GOIDs, and organism, this function will return the species specific gene list annotated by the given GOIDs.
##'
##'
##' @title query genes annotated by given GOIDs
##' @param GOID the query GO IDs
##' @param organism one of human, mouse and yeast.
##' @return a list of gene IDs, the names of the list is the GOIDs
##' @author Guangchuang Yu
getGO2ExtID <- function(GOID, organism) {
    GO2ExtID <- switch(organism,
                       human = mget(GOID, org.Hs.egGO2ALLEGS, ifnotfound=NA),
                       mouse = mget(GOID, org.Mm.egGO2ALLEGS, ifnotfound=NA),
                       yeast = mget(GOID, org.Sc.sgdGO2ALLORFS, ifnotfound=NA),
                       )
    GO2ExtID <- lapply(GO2ExtID, function(i) unique(i))
    return(GO2ExtID)
}

##' query GOIDs at a specific level.
##'
##'
##' @title get GOIDs at a specific level
##' @param ont Ontology
##' @param level GO level
##' @return a vector of GOIDs
##' @author Guangchuang Yu
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

##' generate a bar plot
##'
##' interal use, not for user.
##' @title internal function of barplot
##' @param result a data frame of enrichment result.
##' @param title graph title
##' @param font.size font size
##' @return ggplot object
##' @author Guangchuang Yu
.barplotInternal <- function(result, title, font.size=12) {
    Description <- Count <- NULL # to satisfy codetools
    pg <- ggplot(result, aes(x=Description, y = Count)) +
        geom_bar() +
            coord_flip()
    .pModify(pg, title, font.size)
}

##' changing ggplot object's title and font size
##'
##' interal use, not for user.
##' @title changing title and font size
##' @param p ggplot object
##' @param title graph title
##' @param font.size font size
##' @return ggplot object
##' @author Guangchuang Yu
.pModify <- function(p, title="", font.size=12) {
    p <- p +
        xlab("") +
            ylab("") +
                opts(axis.text.x = theme_text(colour="black", size=font.size, vjust = 1)) +
                    opts(axis.text.y = theme_text(colour="black", size=font.size, hjust = 1)) +
                        opts(title=title)+theme_bw()
    return(p)
}

##' Internal plot function for plotting compareClusterResult
##'
##'
##' @title .PlotClusterProfInternal
##' @param clProf.reshape.df data frame of compareCluster result
##' @param type one of dot and bar
##' @param by one of percentage and count
##' @param title graph title
##' @param font.size graph font size
##' @return ggplot object
##' @author Guangchuang Yu
.PlotClusterProfInternal <- function(clProf.reshape.df,  type = "dot", by = "percentage",title="", font.size=12) {
    Description <- Percentage <- Count <- Cluster <- Pvalue <- pvalue <- NULL # to satisfy codetools
    if (type == "bar") {
        if (by == "percentage") {
            p <- ggplot(clProf.reshape.df, aes(x=Description, y = Percentage, fill=Cluster))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df, aes(x=Description, y = Count, fill=Cluster))
        } else {

        }
        p <- p +
            geom_bar() +
                coord_flip()
    }
    if (type == "dot") {
        if (by == "percentage") {
            p <- ggplot(clProf.reshape.df, aes(x = Cluster, y = Description, size = Percentage))
        } else if (by == "count") {
            p <- ggplot(clProf.reshape.df, aes(x = Cluster, y = Description, size = Count))
        } else {

        }
        if (any(colnames(clProf.reshape.df) == "pvalue")) {
            p <- p +
                geom_point() +
                    aes(color=pvalue) +
                        scale_colour_gradient(low="red", high="yellow")
        } else {
            p <- p + geom_point(colour="steelblue")
        }
    }
    .pModify(p, title, font.size)
}



##' provide a vector of KEGG pathway IDs, this function will convert them to corresponding KEGG pathway Names
##'
##'
##' @title convert KEGG pathway ID to pathway Name
##' @param pathIDs KEGG pathway IDs
##' @return KEGG pathway names
##' @author Guangchuang Yu
path2Name <- function(pathIDs) {
    pathIDs <- gsub("^\\D+", "",pathIDs, perl=T)
    path2name <- mget(pathIDs, KEGGPATHID2NAME)
    return(path2name)
}


##' provide numerator and denominator, return numerator/denominator
##'
##'
##' @title getRatio
##' @param a numerator
##' @param b denominator
##' @return numerator/denominator
##' @author Guangchuang Yu
getRatio <- function(a, b) {
    x=paste(a, "/", b, sep="", collapse="")
    return(x)
}


##' hypergeometric test for enrichment analysis
##'
##'
##' @title hypergeometric test
##' @param numWdrawn number of White balls drawn
##' @param numW number of White balls
##' @param numB number of Black balls
##' @param numDrawn number of balls drawn
##' @return pvalue
##' @author Guangchuang Yu
HyperG <- function(numWdrawn, numW, numB, numDrawn) {
    pvalue <- phyper(numWdrawn, numW, numB, numDrawn, lower.tail=FALSE)
    return(pvalue)
}


##' convert a list of gene IDs to gene Names.
##'
##'
##' @title convert gene IDs to gene Names
##' @param geneID.list a list of gene IDs
##' @param organism one of human, mouse and yeast.
##' @return a list of gene names.
##' @author Guangchuang Yu
geneID2geneName <- function(geneID.list, organism) {
    annotation <- switch(organism,
                         human = org.Hs.egSYMBOL,
                         mouse = org.Mm.egSYMBOL,
                         yeast = org.Sc.sgdGENENAME,
                         )
    gn <- lapply(geneID.list, function(i) unique(unlist(mget(i, annotation))))
    return(gn)
}
