
#' convert goid to descriptive term
#'
#'
#' @title go2term
#' @param goid a vector of GO IDs
#' @return data.frame
#' @export
#' @author Guangchuang Yu
go2term <- function(goid) {
    GO2TERM <- get_GO2TERM_table()
    res <- GO2TERM[GO2TERM[,1] %in% goid, ]
    rownames(res) <- NULL
    return(res)
}

#' convert goid to ontology (BP, CC, MF)
#'
#'
#' @title go2ont
#' @param goid a vector of GO IDs
#' @return data.frame
#' @export
#' @author Guangchuang Yu
go2ont <- function(goid) {
    GO2Ontology <- get_GO2Ontology_table()
    res <- GO2Ontology[GO2Ontology[,1] %in% goid,]
    rownames(res) <- NULL
    return(res)
}


get_GOTERM <- function() {
    pos <- 1
    envir <- as.environment(pos)
    if (!exists(".GOTERM_Env", envir=envir)) {
        assign(".GOTERM_Env", new.env(), envir)
    }
    GOTERM_Env <- get(".GOTERM_Env", envir = envir)
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


excludeGOlevel <- function(x, ont, level) {
    lv <- unlist(lapply(level, getGOLevel, ont=ont))
    x <- excludeGOterm(x, lv)
    return(x)
}

excludeGOterm <- function(x, term) {
    if ( is(x, "enrichResult") ) {
        x@result <- x@result[! x@result[, "ID"] %in% term, ]
    } else if ( is(x, "compareClusterResult") ) {
        x@compareClusterResult <- x@compareClusterResult[! x@compareClusterResult[, "ID"] %in% term, ]
    } else {
        stop("x should be one of enrichResult of compareClusterResult...")
    }
    return(x)
}

keepGOlevel <- function(x, ont, level) {
    lv <- unlist(lapply(level, getGOLevel, ont=ont))
    x <- keepGOterm(x, lv)
    return(x)
}

keepGOterm <- function(x, term) {
    if ( is(x, "enrichResult") ) {
        x@result <- x@result[x@result[, "ID"] %in% term, ]
    } else if ( is(x, "compareClusterResult") ) {
        x@compareClusterResult <- x@compareClusterResult[x@compareClusterResult[, "ID"] %in% term, ]
    } else {
        stop("x should be one of enrichResult of compareClusterResult...")
    }
    return(x)
}

#' query GOIDs at a specific level.
#'
#'
#' @title get GOIDs at a specific level
#' @param ont Ontology
#' @param level GO level
#' @return a vector of GOIDs
#' @importFrom GO.db GOBPCHILDREN
#' @importFrom GO.db GOCCCHILDREN
#' @importFrom GO.db GOMFCHILDREN
#' @importMethodsFrom AnnotationDbi mget
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @noRd
getGOLevel <- function(ont, level) {
    lv <- go_level_map(ont)
    names(lv)[lv %in% level]
}

#' Map every GO term to its level in the ontology
#'
#' Level 1 is the ontology root (`GO:0008150` / `GO:0005575` / `GO:0003674`),
#' level 2 its direct children, and so on, following the `*CHILDREN` graph in
#' GO.db. This is the definition `getGOLevel()` uses to pick the terms at a
#' requested level, and `add_go_level()` uses it to label enrichment results.
#'
#' The returned vector is named by GO ID and may contain an ID more than once:
#' GO is a DAG, so a term can be reachable at several depths and appears in each
#' frontier it belongs to. Callers that want one level per term should keep the
#' first occurrence, which is the shallowest.
#'
#' @param ont Ontology, one of "BP", "CC" or "MF"
#' @param max_level deepest level to walk
#' @return a named integer vector
#' @importFrom GO.db GOBPCHILDREN
#' @importFrom GO.db GOCCCHILDREN
#' @importFrom GO.db GOMFCHILDREN
#' @importFrom stats setNames
#' @importMethodsFrom AnnotationDbi mget
#' @noRd
go_level_map <- function(ont, max_level = 30L) {
    switch(ont,
           MF = {
               topNode <- "GO:0003674"
               Children <- GOMFCHILDREN
           },
           BP = {
               topNode <- "GO:0008150"
               Children <- GOBPCHILDREN
           },
           CC = {
               topNode <- "GO:0005575"
               Children <- GOCCCHILDREN
           },
           stop("ontology should be one of 'MF', 'CC' or 'BP'")
           )

    ids <- topNode
    levels <- 1L
    Node <- topNode
    for (i in seq_len(max_level - 1)) {
        Node <- mget(Node, Children, ifnotfound = NA)
        Node <- unique(as.vector(unlist(Node)))
        Node <- Node[!is.na(Node)]
        if (length(Node) == 0) {
            break
        }
        ids <- c(ids, Node)
        levels <- c(levels, rep(i + 1L, length(Node)))
    }
    setNames(levels, ids)
}

#' Add the GO level of each enriched term
#'
#' `enrichGO()`/`gseGO()` results carry the term, its ontology and its
#' statistics, but not how deep the term sits in the GO hierarchy. This appends
#' that as a `level` column (1 = ontology root), which makes it easy to keep a
#' band of levels by filtering, e.g. `subset(res, level >= 3 & level <= 6)`.
#'
#' The level is taken from the ontology the term itself belongs to, so an
#' `ont = "ALL"` (or `compareCluster()`) result is labelled correctly term by
#' term. Terms that GO.db does not recognise get `NA`.
#'
#' @param x an `enrichResult`, `gseaResult` or `compareClusterResult` from a GO
#'   enrichment analysis
#' @return `x` with a `level` column appended to its result table
#' @export
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @examples
#' \dontrun{
#' x <- enrichGO(gene, OrgDb = "org.Hs.eg.db", ont = "BP")
#' x <- add_go_level(x)
#' subset(x, level >= 3 & level <= 6)
#' }
add_go_level <- function(x) {
    df <- as.data.frame(x)
    if (!"ID" %in% names(df)) {
        stop("the result table has no 'ID' column")
    }
    if (!requireNamespace("GO.db", quietly = TRUE)) {
        stop("GO.db is required to determine GO levels")
    }

    ids <- as.character(df$ID)
    ont_of <- AnnotationDbi::Ontology(GO.db::GOTERM)[ids]
    level <- rep(NA_integer_, length(ids))

    for (ont in unique(ont_of[!is.na(ont_of)])) {
        if (!ont %in% c("BP", "CC", "MF")) {
            next
        }
        ## keep the shallowest level per term: GO is a DAG, so a term can be
        ## reached at more than one depth
        lv <- go_level_map(ont)
        lv <- lv[!duplicated(names(lv))]
        idx <- which(!is.na(ont_of) & ont_of == ont)
        level[idx] <- unname(lv[ids[idx]])
    }

    df$level <- level
    x@result <- df
    return(x)
}


add_GO_Ontology <- function(obj, GO_DATA) {
    if (is(obj, 'gseaResult')) {
        obj@setType <- "GOALL"
    } else if (is(obj, 'enrichResult')) {
        obj@ontology <- 'GOALL'
    }

    df <- obj@result
    
    # Handle GO_DATA as either environment or GSON object
    if (is.environment(GO_DATA)) {
        GO2ONT <- get("GO2ONT", envir=GO_DATA)
    } else if (inherits(GO_DATA, "GSON")) {
        # Extract GO IDs from GSON and get ontology from GO.db
        go_ids <- unique(GO_DATA@gsid2gene$gsid)
        # Get ontology mapping from GO.db
        if (requireNamespace("GO.db", quietly = TRUE)) {
            GO2ONT <- AnnotationDbi::Ontology(GO.db::GOTERM)[go_ids]
            names(GO2ONT) <- go_ids
        } else {
            stop("GO.db package required for ontology mapping")
        }
    } else {
        stop("GO_DATA must be an environment or GSON object")
    }
    
    df <- cbind(ONTOLOGY=GO2ONT[df$ID], df)
    obj@result <- df
    return(obj)
}


get_go_ontology <- function(x) {
    if (is(x, "compareClusterResult")) {
        if (x@fun != "enrichGO" && x@fun != "groupGO" && x@fun != "gseGO") {
            stop("simplify only work for GO...")
        }
        ont <- x@.call$ont
        if (is.null(ont) || !inherits(ont, "character")) {
            ## should be "MF", default value of enrichGO
            ## it's safe to determine from the output
            ont <- x@compareClusterResult$ID[1] %>% GOTERM[[.]] %>% Ontology
        }
    } else if (is(x, "enrichResult")) {
        if (!x@ontology %in% c("BP", "MF", "CC"))
            stop("ontology should be one of 'MF', 'BP', 'CC'...")

        ont <- x@ontology
    } else {
        stop("x should be an instance of 'enrichResult' or 'compareClusterResult'...")
    }

    return(ont)
}
