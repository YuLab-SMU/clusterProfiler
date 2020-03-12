
##' @method arrange enrichResult
##' @export
arrange.enrichResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% arrange(!!!dots,)
    return(.data)
}

##' @method arrange gseaResult
##' @export
arrange.gseaResult <- arrange.enrichResult


##' @method arrange compareClusterResult
##' @export
arrange.compareClusterResult <- function(.data, ...) {
    dots <- quos(...)
    .data@compareClusterResult %<>% arrange(!!!dots,)
    return(.data)
}
