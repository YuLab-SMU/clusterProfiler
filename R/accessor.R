##' @method [ compareClusterResult
##' @export
`[.compareClusterResult` <- function(x, i, j) {
              x@result[i,j]
}

##' @method [[ compareClusterResult
##' @export
`[[.compareClusterResult` <- function(x, i) {
    if (!i %in% names(x@geneInCategory))
        stop("input term not found...")
    x@geneInCategory[[i]]
}


##' @importFrom utils head
##' @method head compareClusterResult
##' @export
head.compareClusterResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}


##' @importFrom utils tail
##' @method tail compareClusterResult
##' @export
tail.compareClusterResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}


##' @method dim compareClusterResult
##' @export
dim.compareClusterResult <- function(x) {
    dim(x@result)
}
