dotplot.compareClusterResult <- function(object, x=~Cluster, colorBy="p.adjust", showCategory=5, by="geneRatio", includeAll=TRUE, font.size=12, title="") {
    df <- fortify(object, showCategory=showCategory, by=by, includeAll=includeAll)
    plotting.clusterProfile(df, x=x, type="dot", colorBy=colorBy, by=by, title=title, font.size=font.size)
}

##' convert compareClusterResult to a data.frame that ready for plot
##'
##'
##' @rdname fortify
##' @title fortify
##' @param model compareClusterResult object
##' @param data not use here
##' @param showCategory category numbers
##' @param by one of geneRatio, Percentage or count
##' @param includeAll logical
##' @return data.frame
##' @importFrom ggplot2 fortify
##' @importFrom plyr ddply
##' @importFrom plyr mdply
##' @importFrom plyr .
##' @method fortify compareClusterResult
##' @export
##' @author Guangchuang Yu
fortify.compareClusterResult <- function(model, data, showCategory=5, by="geneRatio", includeAll=TRUE) {
    clProf.df <- as.data.frame(model)

    ## get top 5 (default) categories of each gene cluster.
    if (is.null(showCategory)) {
        result <- clProf.df
    } else {
        Cluster <- NULL # to satisfy codetools
        result <- ddply(.data = clProf.df,
                        .variables = .(Cluster),
                        .fun = function(df, N) {
                            if (length(df$Count) > N) {
                                if (any(colnames(df) == "pvalue")) {
                                    idx <- order(df$pvalue, decreasing=FALSE)[1:N]
                                } else {
                                    ## for groupGO
                                    idx <- order(df$Count, decreasing=T)[1:N]
                                }
                                return(df[idx,])
                            } else {
                                return(df)
                            }
                        },
                        N=showCategory
                        )

    }
    ID <- NULL
    if (includeAll == TRUE) {
        result = subset(clProf.df, ID %in% result$ID)
    }

    ## remove zero count
    result$Description <- as.character(result$Description) ## un-factor
    GOlevel <- result[,c("ID", "Description")] ## GO ID and Term
    GOlevel <- unique(GOlevel)

    result <- result[result$Count != 0, ]
    result$Description <- factor(result$Description,
                                 levels=rev(GOlevel[,2]))


    if (by=="rowPercentage") {
        Description <- Count <- NULL # to satisfy codetools
        result <- ddply(result,
                        .(Description),
                        transform,
                        Percentage = Count/sum(Count),
                        Total = sum(Count))

        ## label GO Description with gene counts.
        x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
        y <- sapply(x[,3], paste, ")", sep="")
        result$Description <- y

        ## restore the original order of GO Description
        xx <- result[,c(2,3)]
        xx <- unique(xx)
        rownames(xx) <- xx[,1]
        Termlevel <- xx[as.character(GOlevel[,1]),2]

        ##drop the *Total* column
        result <- result[, colnames(result) != "Total"]

        result$Description <- factor(result$Description,
                                     levels=rev(Termlevel))

    } else if (by == "count") {
        ## nothing
    } else if (by == "geneRatio") {
        gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
        gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
        result$GeneRatio = gsize/gcsize
        result$Cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")
    } else {
        ## nothing
    }
    return(result)
}
