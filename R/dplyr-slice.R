##' @method slice enrichResult
##' @export 
slice.enrichResult <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@result %<>% slice(!!!dots, .preserve = .preserve)
    return(.data)
}

##' @method slice gseaResult
##' @export
slice.gseaResult <- slice.enrichResult

##' @method slice compareClusterResult
##' @export 
slice.compareClusterResult <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@compareClusterResult %<>% slice(!!!dots, .preserve = .preserve)
    return(.data)
}

