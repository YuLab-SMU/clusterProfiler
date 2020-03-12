##' @method mutate enrichResult
##' @importFrom dplyr mutate
##' @importFrom rlang quos
##' @export
##' @author Guangchuang Yu 
mutate.enrichResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% mutate(!!!dots)
    return(.data)
}

##' @method mutate gseaResult
##' @export
mutate.gseaResult <- mutate.enrichResult


##' @method mutate compareClusterResult
##' @export
mutate.compareClusterResult <- function(.data, ...) {
    dots <- quos(...)
    .data@compareClusterResult %<>% mutate(!!!dots)
    return(.data)
}
