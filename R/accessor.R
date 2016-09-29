##' @method as.data.frame compareClusterResult
##' @export
as.data.frame.compareClusterResult <- function(x, ...) {
    as.data.frame(x@compareClusterResult, ...)
}

##' @method [ compareClusterResult
##' @export
`[.compareClusterResult` <- function(x, i, j) {
              x@result[i,j]
}

##' @method [[ compareClusterResult
##' @export
`[[.compareClusterResult` <- function(x, i) {
    gc <- geneInCategory(x)
    if (!i %in% names(gc))
        stop("input term not found...")
    gc[[i]]
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


##' @method geneID groupGOResult
##' @importFrom DOSE geneID
##' @export
geneID.groupGOResult <- function(x) as.character(x@result$geneID)


##' @method geneInCategory groupGOResult
##' @export
##' @importFrom DOSE geneInCategory
##' @importFrom stats setNames
geneInCategory.groupGOResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))
