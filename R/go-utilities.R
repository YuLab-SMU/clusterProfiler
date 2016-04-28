
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

##' convert goid to ontology (BP, CC, MF)
##'
##' 
##' @title go2ont
##' @param goid a vector of GO IDs
##' @return data.frame
##' @export
##' @author Guangchuang Yu
go2ont <- function(goid) {
    GO2Ontology <- get_GO2Ontology_table()
    res <- GO2Ontology[GO2Ontology[,1] %in% goid,]
    rownames(res) <- NULL
    return(res)
}


get_GOTERM <- function() {
    if (!exists(".GOTERM_Env", envir=.GlobalEnv)) {
        assign(".GOTERM_Env", new.env(), .GlobalEnv)
    }
    GOTERM_Env <- get(".GOTERM_Env", envir = .GlobalEnv)
    if (exists("GOTERM.df", envir = GOTERM_Env)) {
        GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
    } else {
        GOTERM.df <- toTable(GOTERM)
        assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
    }
    return(GOTERM.df)
}

get_GO2TERM_table <- function() {
    GOTERM.df <- get_GOTERM()
    GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GO2Ontology_table <- function() {
    GOTERM.df <- get_GOTERM()
    GOTERM.df[, c("go_id", "Ontology")] %>% unique
}
