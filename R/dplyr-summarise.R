##' @method summarise enrichResult
##' @export
summarise.enrichResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %>% summarise(!!!dots)
}


##' @method summarise gseaResult
##' @export
summarise.gseaResult <- summarise.enrichResult

##' @method summarise compareClusterResult
##' @export
summarise.compareClusterResult <- function(.data, ...) {
    dots <- quos(...)
    .data@compareClusterResult %>% summarise(!!!dots)
}
