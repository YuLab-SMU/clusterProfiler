##' @method filter enrichResult
##' @export
filter.enrichResult <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@result %<>% filter(!!!dots, .preserve = .preserve)
    return(.data)
}


##' @method filter gseaResult
##' @export
filter.gseaResult <- filter.enrichResult

##' @method filter compareClusterResult
##' @export
filter.compareClusterResult <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@compareClusterResult %<>% filter(!!!dots, .preserve = .preserve)
    return(.data)
}

