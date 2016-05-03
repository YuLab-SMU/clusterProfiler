get_kegg_species <- function() {
    x <- XML::readHTMLTable("http://www.genome.jp/kegg/catalog/org_list.html") 
    y <- get_species_name(x[[2]], "Eukaryotes")
    y2 <- get_species_name(x[[3]], 'Prokaryotes')
    
    sci_name <- gsub(" \\(.*$", '', y[,2])
    com_name <- gsub("[^\\(]+ \\(([^\\)]+)\\)$", '\\1', y[,2])
    eu <- data.frame(kegg_code=unlist(y[,1]),
                     scientific_name = sci_name,
                     common_name = com_name, 
                     stringsAsFactors = FALSE)
    pr <- data.frame(kegg_code=unlist(y2[,1]),
                     scientific_name = unlist(y2[,2]),
                     common_name = NA,
                     stringsAsFactors = FALSE)
    kegg_species <- rbind(eu, pr)
    save(kegg_species, file="kegg_species.rda")
    invisible(kegg_species)
}

get_species_name <- function(y, table) {
    idx <- get_species_name_idx(y, table)
    t(sapply(1:nrow(idx), function(i) {
        y[] = lapply(y, as.character)
        y[i, idx[i,]]
    }))
}


get_species_name_idx <- function(y, table='Eukaryotes') {
    table <- match.arg(table, c("Eukaryotes", "Prokaryotes"))
    t(apply(y, 1, function(x) {
        ii <- which(!is.na(x))
        n <- length(ii)
        if (table == "Eukaryotes") {
            return(ii[(n-2):(n-1)])
        } else {
            return(ii[(n-3):(n-2)])
        }
    }))
}


