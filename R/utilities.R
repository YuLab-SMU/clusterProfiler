.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
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
##' @return ggplot object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 %+%
##' @importFrom ggplot2 opts
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme_text
##' @importFrom ggplot2 scale_colour_gradient
##' @author Guangchuang Yu \url{http://ygc.name}
plotting.clusterProfile <- function(clProf.reshape.df,  type = "dot", by = "percentage",title="", font.size=12) {
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
                        scale_colour_gradient(low="red", high="blue")
        } else {
            p <- p + geom_point(colour="steelblue")
        }
    }
    p <- p + xlab("") + ylab("") +
        opts(axis.text.x = theme_text(colour="black", size=font.size, vjust = 1)) +
            opts(axis.text.y = theme_text(colour="black", size=font.size, hjust = 1)) +
                opts(title=title)+theme_bw()
    return(p)
}

