.initial <- function() {
    assign("clusterProfilesEnv", new.env(),.GlobalEnv)
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
                       ecolik12     = "org.EcK12.egGO2ALLEGS",
                       ecsakai      = "org.EcSakai.egGO2ALLEGS",
                       fly          = "org.Dm.egGO2ALLEGS",
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
                                    font.size=12,
                                    angle.axis.x=90,
                                    hjust.axis.x=1,
                                    vjust.axis.x=0.5) {
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
    p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
                       hjust=hjust.axis.x,
                       vjust=vjust.axis.x))
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
