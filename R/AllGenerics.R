##' @importFrom BiocGenerics plot
##' @exportMethod plot
if ( !isGeneric("plot") )
	setGeneric("plot", function(x, ...) standardGeneric("plot"))
