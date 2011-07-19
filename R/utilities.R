
##' Assign some initial instances.
##'
##'
##' @title initial function for clusterProfiler
##' @return environment instance
##' @author Guangchuang Yu
.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param GOID
##' @return GO Terms
##' @author Guangchuang Yu
GO2Term <- function(GOID) {
    go <- mget(GOID, GOTERM, ifnotfound=NA)
    term <- sapply(go, Term)
    return(term)
}

getGO2ExtID <- function(Terms, organism) {
    GO2ExtID <- switch(organism,
                       human = mget(Terms, org.Hs.egGO2ALLEGS, ifnotfound=NA),
                       mouse = mget(Terms, org.Mm.egGO2ALLEGS, ifnotfound=NA),
                       yeast = mget(Terms, org.Sc.sgdGO2ALLORFS, ifnotfound=NA),
                       )
    GO2ExtID <- lapply(GO2ExtID, function(i) unique(i))
    return(GO2ExtID)
}


getGOLevel <- function(ont, level) {
    ## get GO nodes at a specific level...
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

.barplotInternal <- function(result, title, font.size=12) {
    Description <- Count <- NULL # to satisfy codetools
    pg <- ggplot(result, aes(x=Description, y = Count)) +
        geom_bar() +
            coord_flip()
    .pModify(pg, title, font.size)
}

.pModify <- function(p, title="", font.size=12) {
    p <- p + xlab("") + ylab("") + opts(axis.text.x = theme_text(colour="black", size=font.size, vjust = 1)) + opts(axis.text.y = theme_text(colour="black", size=font.size, hjust = 1)) + opts(title=title)+theme_bw()
    return(p)
}

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






#### KEGG Enrichment Analysis ###
path2Name <- function(pathIDs) {
    pathIDs <- gsub("^\\D+", "",pathIDs, perl=T)
    path2name <- mget(pathIDs, KEGGPATHID2NAME)
    return(path2name)
}

getRatio <- function(a, b) {
    x=paste(a, "/", b, sep="", collapse="")
    return(x)
}

HyperG <- function(numWdrawn, numW, numB, numDrawn) {
    ##numWdrawn: number of White balls drawn
    ##numW: number of White balls
    ##numB: number of Black balls
    ##numDrawn: number of balls drawn
    pvalue <- phyper(numWdrawn, numW, numB, numDrawn, lower.tail=FALSE)
    return(pvalue)
}

geneID2geneName <- function(geneID.list, organism) {
    annotation <- switch(organism,
                         human = org.Hs.egSYMBOL,
                         mouse = org.Mm.egSYMBOL,
                         yeast = org.Sc.sgdGENENAME,
                         )
    gn <- lapply(geneID.list, function(i) unique(unlist(mget(i, annotation))))
    return(gn)
}
