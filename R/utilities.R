
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
##' @importFrom ggplot2 scale_color_continuous
##' @importFrom ggplot2 guide_colorbar
##' @importFrom DOSE theme_dose
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
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
                scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
            ## scale_color_gradientn(guide=guide_colorbar(reverse=TRUE), colors = enrichplot:::sig_palette)
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


GSEA_internal <- DOSE:::GSEA_internal
enricher_internal <- DOSE:::enricher_internal
get_fun_from_pkg <- rvcheck:::get_fun_from_pkg

