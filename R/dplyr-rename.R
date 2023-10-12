
##' @method rename enrichResult
##' @export
rename.enrichResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% rename(!!!dots,)
    return(.data)
}

##' @method rename gseaResult
##' @export
rename.gseaResult <- rename.enrichResult

##' @method rename compareClusterResult
##' @export
rename.compareClusterResult <- function(.data, ...) {
    dots <- quos(...)
    .data@compareClusterResult %<>% rename(!!!dots,)
    return(.data)
}
