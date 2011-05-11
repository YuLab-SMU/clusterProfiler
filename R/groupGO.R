groupGO <- function(gene, organism="human", ont="CC", level = 2) { 
	GOLevel <- .getGOLevel(ont, level)
	goGene= .goGene(GOLevel, gene, organism)
	GeneIDs = .getGeneID(goGene)
	Count <- unlist(lapply(goGene, length))
	Descriptions <- sapply(GOLevel, .GO2Term)
	result = data.frame(GOID=GOLevel, Description=Descriptions, Count=Count, GeneID=GeneIDs)
	new("groupGOResult", 
		groupGOResult=result, 
		Ont = ont,
		Level = level,
		Organism = organism,
		Gene = gene
	)
}

setClass("groupGOResult",
	representation=representation(
		groupGOResult="data.frame",
		Ont = "character",
		Level = "numeric",
		Organism = "character",
		Gene = "character"
	)
)
	
setMethod("show", signature(object="groupGOResult"),
	function (object){
		ont = object@Ont
		Level = object@Level
		Organism = object@Organism
		Gene = object@Gene
		cat ("GO", ont, "Profiles", "at level", Level, "of", length(Gene),  Organism, "genes", "\n")
	}
)

setMethod("summary", signature(object="groupGOResult"),
	function (object){
		return(object@groupGOResult)
	}
)

setMethod("plot", signature(x="groupGOResult"),
	function (x, order="FALSE", caption="", font.size=12, drop=FALSE){
		groupGOResult <- summary(x)
		if (drop == TRUE) {
			# drop void category.
			groupGOResult <- groupGOResult[groupGOResult$Count != 0, ]
		}
		if (order == TRUE) {
			idx <- order(groupGOResult$Count)
			groupGOResult <- groupGOResult[idx,]
		}
		groupGOResult$Description <- factor(groupGOResult$Description, level= as.character(groupGOResult$Description))
		.barplotInternal(groupGOResult, caption, font.size)
	}
)