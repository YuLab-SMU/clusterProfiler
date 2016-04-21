
##' convert goid to descriptive term
##'
##' 
##' @title go2term
##' @param goid a vector of GO IDs
##' @return data.frame
##' @export
##' @author Guangchuang Yu
go2term <- function(goid) {
    GO2TERM <- get_GO2TERM_table()
    res <- GO2TERM[GO2TERM[,1] %in% goid, ]
    rownames(res) <- NULL
    return(res)
}

get_GO2TERM_table <- function() {
    if (!exists(".GO2TERM_Env", envir=.GlobalEnv)) {
        assign(".GO2TERM_Env", new.env(), .GlobalEnv)
    }
    GO2TERM_Env <- get(".GO2TERM_Env", envir = .GlobalEnv)
    if (exists("GO2TERM", envir = GO2TERM_Env)) {
        GO2TERM <- get("GO2TERM", envir=GO2TERM_Env)
    } else {
        goids <- toTable(GOTERM)
        GO2TERM <- goids[, c("go_id", "Term")] %>% unique
        assign("GO2TERM", GO2TERM, envir = GO2TERM_Env)
    }
    return(GO2TERM)
}
