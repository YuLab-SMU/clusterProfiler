
##' @method group_by enrichResult
##' @export
group_by.enrichResult <- function(.data, ..., add = FALSE, .drop = FALSE) {
    dots <- quos(...)
    .data@result %<>% group_by(!!!dots, add = add, .drop = .drop)
    return(.data)
}


##' @method group_by gseaResult
##' @export
group_by.gseaResult <- group_by.enrichResult

##' @method group_by compareClusterResult
##' @export
group_by.compareClusterResult <- function(.data, ..., add = FALSE, .drop = FALSE) {
    dots <- quos(...)
    .data@compareClusterResult %<>% group_by(!!!dots, add = add, .drop = .drop)
    return(.data)
}
