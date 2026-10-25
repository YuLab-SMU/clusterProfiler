#' Record evidence-source failures instead of swallowing them
#'
#' @description
#' Every source this package consults can fail: a database rejects an
#' oversized request, an endpoint is unreachable, a companion package is not
#' installed. Those failures used to be caught with `error = function(e) NULL`,
#' which made them indistinguishable from a successful query that found nothing.
#'
#' The distinction matters because the two lead to opposite conclusions. "No
#' knockdown dataset moves these genes" is a finding. "GPSAdb rejected the
#' request" is not, and reporting the second as the first produces a confident
#' negative that nothing in the output contradicts.
#'
#' The rule this file implements is:
#'
#' * A source that fails becomes **`failed`**, which is neither `absent` nor
#'   `untested`, and is reported by `print()` and in `$diagnostics$failures`.
#' * If **every** consulted source fails, the call aborts: there is nothing to
#'   report a verdict about.
#' * A **check** that cannot run aborts rather than passing. A boundary that
#'   silently admits everything is worse than no boundary, because the caller
#'   believes it was enforced.
#'
#' @name mechanism-failures
#' @keywords internal
NULL

## One failure. Kept small and flat so it can be row-bound cheaply and printed
## as a table.
mechanism_failure <- function(source, stage, message, fatal = FALSE) {
    data.frame(
        source  = as.character(source),
        stage   = as.character(stage),
        message = as.character(message),
        fatal   = isTRUE(fatal),
        stringsAsFactors = FALSE
    )
}

.mechanism_failures <- new.env(parent = emptyenv())

## Collect failures during one call, then hand them to the bundle.
##
## A collector rather than a global: two concurrent evidence() calls must not
## see each other's failures.
mechanism_failure_collector <- function() {
    e <- new.env(parent = emptyenv())
    e$rows <- list()
    e$add <- function(...) {
        e$rows[[length(e$rows) + 1L]] <- mechanism_failure(...)
        invisible(NULL)
    }
    e$all <- function() {
        if (!length(e$rows)) {
            ## Build the empty frame directly: mechanism_failure() is scalar, and
            ## calling it with zero-length arguments makes data.frame() complain
            ## about mismatched rows.
            return(structure(
                data.frame(source = character(), stage = character(),
                           message = character(), fatal = logical(),
                           stringsAsFactors = FALSE),
                class = c("mechanismFailures", "data.frame")))
        }
        out <- do.call(rbind, e$rows)
        structure(out, class = c("mechanismFailures", "data.frame"))
    }
    e$n <- function() length(e$rows)
    e
}

## Run `expr`, recording a failure instead of swallowing it.
##
## `on_fail` is what the caller wants back so the computation can continue; it
## is never the answer, only a placeholder, and every call site that uses it is
## expected to surface the failure.
mechanism_try <- function(collector, source, stage, expr, on_fail = NULL) {
    tryCatch(
        expr,
        error = function(e) {
            if (!is.null(collector)) {
                collector$add(source = source, stage = stage,
                              message = conditionMessage(e))
            }
            on_fail
        }
    )
}

## Abort when nothing worked.
##
## A bundle whose every source failed would otherwise report an absence for
## every entity, which reads as "the data says no" when the truth is "the data
## was never obtained".
mechanism_abort_if_all_failed <- function(collector, sources) {
    if (is.null(collector) || !collector$n()) return(invisible(TRUE))
    rows <- collector$all()
    failed <- unique(rows$source)
    if (!all(sources %in% failed)) return(invisible(TRUE))

    detail <- paste0("  * ", rows$source, " [", rows$stage, "]: ", rows$message,
                     collapse = "\n")
    rlang::abort(c(
        "Every evidence source failed, so there is nothing to report a verdict about.",
        x = paste0("Failed: ", paste(failed, collapse = ", ")),
        i = paste0("Details:\n", detail),
        i = "Returning an empty bundle would read as 'the evidence says no'."
    ))
}

#' @rdname mechanism-failures
#' @param x An object carrying failures.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.mechanismFailures <- function(x, ...) {
    if (!nrow(x)) {
        cat("<mechanismFailures: none>\n")
        return(invisible(x))
    }
    cat("<mechanismFailures: ", nrow(x), ">\n", sep = "")
    for (i in seq_len(nrow(x))) {
        cat("  ", x$source[i], " [", x$stage[i], "]: ", x$message[i], "\n", sep = "")
    }
    invisible(x)
}
