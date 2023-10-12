##' @method select enrichResult
##' @importFrom dplyr select
##' @export
select.enrichResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% select(!!!dots,)
    return(.data)
}

##' @method select gseaResult
##' @export
select.gseaResult <- select.enrichResult


##' @method select compareClusterResult
##' @export
select.compareClusterResult <- function(.data, ...) {
    dots <- quos(...)
    .data@compareClusterResult %<>% select(!!!dots,)
    return(.data)
}
