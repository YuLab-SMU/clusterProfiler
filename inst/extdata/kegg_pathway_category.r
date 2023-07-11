library(rvest)

url <- 'https://www.genome.jp/kegg/pathway.html'
webpage <- read_html(url)
categories <- html_nodes(webpage, "b") |> html_text()
#categories <- categories[grep("^\\d", categories)]
cat1 <- categories[grep("^\\d\\.\\s", categories)]
cat2 <- categories[grep("^\\d\\.\\d", categories)]

cat1 <- cat1[as.numeric(sub("(\\d).*", "\\1", cat2))]
cat1 <- sub("^\\d\\.\\s", "", cat1)
cat2 <- sub("^\\d\\.\\d+\\s", "", cat2)

keggmap <- html_nodes(webpage, "dl") |> html_text()
keggmap <- keggmap[grep("^\\d{5}", keggmap)]

term <- lapply(keggmap, function(x) {
    y <- unlist(strsplit(x, split="\n"))
    y <- sub("^\\s*", "", y)
    y <- y[y!= ""]

    ## remove these indicators
    ##
    ## KEGG PATHWAY is integrated with MODULE and NETWORK databases as indicated below.
    ## M - module
    ## R - reaction module
    ## N - network
    ##

    y <- sub("^(\\d+)\\sM\\sN(.*)$", "\\1\\2", y)
    y <- sub("^(\\d+)\\sM\\sR(.*)$", "\\1\\2", y)
    y <- sub("^(\\d+)\\sM(.*)$", "\\1\\2", y)
    y <- sub("^(\\d+)\\sR(.*)$", "\\1\\2", y)
    y <- sub("^(\\d+)\\sN(.*)$", "\\1\\2", y)

    id <- sub("^(\\d+)\\D.*$", "\\1", y)
    name <- sub("^\\d+", "", y)
    d <- data.frame(id = id, name = name)

    ii <- grep("Including:", name)
    if (length(ii) > 0) {
        for (i in ii) {
            tmp <- unlist(strsplit(name[i], split="Including: "))
            d$name[i] <- tmp[1]
            d2 <- data.frame(
                id = id[i],
                name = unlist(strsplit(tmp[2], split=", "))
            )
            d <- rbind(d, d2)
        }    
    }
    return(d)
})


cat.df <- tibble::tibble(
    category = cat1,
    subcategory = cat2,
    term = term
)
kegg_category <- tidyr::unnest(cat.df, cols = c("term"))

cat(sprintf("--> Number of KEGG pathways: %s\n", nrow(kegg_category)))

save(kegg_category, file="data/kegg_category.rda")
