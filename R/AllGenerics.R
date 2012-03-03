##' plot method generics
##'
##'
##' @docType methods
##' @name plot
##' @rdname plot-methods
##' @title plot method
##' @param ... Additional argument list
##' @return plot
##' @importFrom graphics plot
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
if ( !isGeneric("plot") )
	setGeneric("plot", function(x, ...) standardGeneric("plot"))
