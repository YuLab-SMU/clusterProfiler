########### utilities.R
.initial <- function() {
	assign("clusterProfilesEnv", new.env(),.GlobalEnv)
}

.GO2Term <- function(GOID) {
	go <- get(GOID, GOTERM)
	term <- go@Term
	return(term)
}

.goGene <- function(GOID, gene, organism="human") {
	# return a list of *GO* annotated in *gene* list.
	#######
	if (organism == "human") {
		goAllGene= mget(GOID, org.Hs.egGO2ALLEGS, ifnotfound=NA)
	} else if (organism == "mouse") {
		goAllGene= mget(GOID, org.Mm.egGO2ALLEGS, ifnotfound=NA)
	} else if (organism == "yeast") {
		goAllGene= mget(GOID, org.Sc.sgdGO2ALLORFS, ifnotfound=NA)
	} else {
		stop("only *human* , *mouse*, and *yeast* supported right now...")
	}
	#########
	goGene <- lapply(goAllGene, function(x) gene[gene %in% x])
	return(goGene)
}

.getGeneID <- function(goGene) {
	### return a vector of GeneID list from the output of function *.goGene*
	GeneIDs <- lapply(goGene, function(x) paste(x, collapse="/"))
	GeneIDs <- unlist(GeneIDs)
	return(GeneIDs)
}

.getGOLevel <- function(ont, level) {
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

.barplotInternal <- function(result, caption, font.size=12) {
	Description <- Count <- NULL # to satisfy codetools
	pg <- ggplot(result, aes(x=Description, y = Count)) + geom_bar() + coord_flip()
	.pModify(pg, caption, font.size)
}

.pModify <- function(p, caption="", font.size=12) {
	p <- p + xlab("") + ylab("") + opts(axis.text.x = theme_text(colour="black", size=font.size, vjust = 1)) + opts(axis.text.y = theme_text(colour="black", size=font.size, hjust = 1)) + opts(title=caption)+theme_bw()
	return(p)
}

.PlotClusterProfInternal <- function(clProf.reshape.df,  type = "dot", by = "percentage",caption="", font.size=12) {
	Description <- Percentage <- Count <- Cluster <- Pvalue <- pvalue <- NULL # to satisfy codetools
	if (type == "bar") {
		if (by == "percentage") {
			p <- ggplot(clProf.reshape.df, aes(x=Description, y = Percentage, fill=Cluster))
		} else if (by == "count") {
			p <- ggplot(clProf.reshape.df, aes(x=Description, y = Count, fill=Cluster))
		} else {

		}
		p <- p+ geom_bar() + coord_flip()
	}
	if (type == "dot") {
		if (by == "percentage") {
			p <- ggplot(clProf.reshape.df, aes(x = Cluster, y = Description, size = Percentage))
		} else if (by == "count") {
			p <- ggplot(clProf.reshape.df, aes(x = Cluster, y = Description, size = Count))
		} else {

		}
		if (any(colnames(clProf.reshape.df) == "Pvalue")) {
			p <- p + geom_point() + aes(color=Pvalue)
			# + scale_colour_gradient(low="red", high="yellow")
		} else if (any(colnames(clProf.reshape.df) == "pvalue")) {
			p <- p + geom_point() + aes(color=pvalue)
		} else {
			p <- p + geom_point(colour="steelblue")
		}
	}
	.pModify(p, caption, font.size)
}


.compareClusterResultTopN <- function(clusterProfResult, limit) {
	Cluster <- NULL # to satisfy codetools
	if (is.null(limit)) {
		return(summary(clusterProfResult))
	}
	clProf.df <- summary(clusterProfResult)
	ddply(clProf.df, .(Cluster), .topN, N=limit)
}

.topN <- function(df, N) {
	if (length(df$Count) > N) {
		idx <- order(df$Count, decreasing=T)[1:N]
		return(df[idx,])
	} else {
		return(df)
	}
}

.compareClusterResultReshape <- function(clProf.df) {
	Description <- Count <- NULL # to satisfy codetools
	GOlevel <- clProf.df[,c(2,3)]
	GOlevel <- unique(GOlevel)

	clProf.df <- ddply(clProf.df, .(Description), transform, Percentage = Count/sum(Count), Total = sum(Count))
	#clProf.df <- ddply(clProf.df, .(Description), transform, Total = sum(Count))

	x <- mdply(clProf.df[, c("Description", "Total")], paste, sep=" (")
	y <- sapply(x[,3], paste, ")", sep="")
	clProf.df$Description <- y		### label GO Description with gene counts.
	xx <- clProf.df[,c(2,3)]
	xx <- unique(xx)
	rownames(xx) <- xx[,1]
	Termlevel <- xx[as.character(GOlevel[,1]),2]

	#clProf.df <-  clProf.df[, -6] ###drop the *Total* column##

	clProf.df$Description <- factor(clProf.df$Description, levels=rev(Termlevel))
	return(clProf.df)
}

.removeZeroCount <- function(clProf.df) {
	## un-factor
	clProf.df$Description <- as.character(clProf.df$Description)
	GOlevel <- clProf.df[,c(2,3)]
	GOlevel <- unique(GOlevel)

	clProf.df <- clProf.df[clProf.df$Count != 0, ]
	clProf.df$Description <- factor(clProf.df$Description, levels=rev(GOlevel[,2]))
	return(clProf.df)
}


#### KEGG Enrichment Analysis ###
path2Name <- function(pathIDs) {
	pathIDs <- gsub("^\\D+", "",pathIDs, perl=T)
	path2name <- mget(pathIDs, KEGGPATHID2NAME)
	return(path2name)
}

.yPaste <- function(a, b) {
	x=paste(a, "/", b, sep="", collapse="")
	return(x)
}

HyperG <- function(numWdrawn, numW, numB, numDrawn) {
	#numWdrawn: number of White balls drawn
	#numW: number of White balls
	#numB: number of Black balls
	#numDrawn: number of balls drawn
	pvalue <- phyper(numWdrawn, numW, numB, numDrawn, lower.tail=FALSE)
	return(pvalue)
}

geneID2geneName <- function(geneID.list, organism) {
    annotation <- switch(organism,
                         human = org.Hs.egSYMBOL,
                         mouse = org.Mm.egSYMBOL,
                         yeast = org.Sc.sgdGENENAME,
                         )
    gn <- lapply(geneID.list, function(i) unique(unlist(mget(i, org.Hs.egSYMBOL))))
    return(gn)
}
