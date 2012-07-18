.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
}



getSupported_Org <- function() {
    supported_Org <- c("human", "mouse", "yeast", "zebrafish")
    return(supported_Org)
}

getAnnoDb <- function(organism) {
    annoDb <- switch(organism,
                     human = "org.Hs.eg.db",
                     mouse = "org.Mm.eg.db",
                     yeast = "org.Sc.sgd.db",
                     zebrafish = "org.Dr.eg.db"
                     )
    return(annoDb)
}

getGO2ALLEG_MappedDb <- function(organism) {
    annoDb <- getAnnoDb(organism)
    require(annoDb, character.only = TRUE)

    mappedDb <- switch(organism,
                       human = "org.Hs.egGO2ALLEGS",
                       mouse = "org.Mm.egGO2ALLEGS",
                       yeast = "org.Sc.sgdGO2ALLORFS",
                       zebrafish = "org.Dr.egGO2ALLEGS"
                       )
    mappedDb <- eval(parse(text=mappedDb))
    return(mappedDb)
}

getEG2GO_MappedDb <- function(organism) {
    annoDb <- getAnnoDb(organism)

    require(annoDb, character.only = TRUE)

    mappedDb <- switch(organism,
                       human = "org.Hs.egGO",
                       mouse = "org.Mm.egGO",
                       yeast = "org.Sc.sgdGO",
                       zebrafish = "org.Dr.egGO"
                       )
    mappedDb <- eval(parse(text=mappedDb))
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

##' building GO mapping files
##'
##' provided by a data.frame of gene and GO directly annotation file
##' this function will building gene to GO and GO to gene mapping,
##' with directly and undirectly annotation.
##' @title buildGOmap
##' @param gomap data.frame with two columns names "entrezgene", and "go_accession"
##' @return files save in the the working directory
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
##' @importFrom plyr dlply
##' @export
##' @author Yu Guangchuang
buildGOmap <- function(gomap) {
    if( any( colnames(gomap) %in% "go_id" ) ) {
        colnames(gomap)[colnames(gomap) %in% "go_id"] <- "go_accession"
    }

    ## remove empty GO annotation
    gomap <- gomap[gomap$go_accession != "", ]


    GO2EG <- dlply(gomap,"go_accession",.fun=function(i) i$entrezgene)
    EG2GO <- dlply(gomap,"entrezgene",.fun=function(i) i$go_accession)

    save(GO2EG, file="GO2EG.rda", compress="xz")
    save(EG2GO, file="EG2GO.rda", compress="xz")



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
    save(EG2ALLGO, file="EG2ALLGO.rda", compress="xz")

    len <- lapply(EG2ALLGO,length)
    EG2ALLGO.df <- data.frame(EG=rep(names(EG2ALLGO), times=len),
                              GO=unlist(EG2ALLGO))
    GO <- NULL ## satisfy code tools
    GO2ALLEG <- dlply(EG2ALLGO.df, .(GO), function(i) as.character(i$EG))
    GO2ALLEG <- lapply(GO2ALLEG, unique)
    save(GO2ALLEG, file="GO2ALLEG.rda", compress="xz")
    print("GO Annotation Mapping files save in the working directory.")
}
