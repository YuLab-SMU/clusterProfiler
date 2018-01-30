##' retreve annotation data from uniprot
##'
##'
##' @title uniprot_get
##' @param taxID taxonomy ID
##' @return gene table data frame
##' @export
##' @author guangchuang yu
uniprot_get <- function(taxID) {
    url <- paste("http://www.uniprot.org/uniprot/?query=taxonomy%3a", taxID,
                 "&force=yes&format=txt", sep="", collapse="")
    parseEmbl(url)
}

##' @importFrom plyr ldply
parseEmbl <- function(file) {
    x <- readLines(file)
    idx <- grep("^//", x)
    if (idx[length(idx)] == length(x)) {
        idx <- idx[-length(idx)]
    }
    begin <- c(1, idx+1)
    end <- c(idx, length(x)) - 1

    res <- lapply(1:length(begin), function(i) {
        item <- x[begin[i]:end[i]]
        parseEmblItem(item)
    })
    res <- ldply(res)
    return(res)
}

parseEmblItem <- function(item) {
    ## uniprot accession
    ac <- item[grep("^AC", item)]
    ac <- gsub("^AC\\s+", "", ac)
    ac <- gsub("\\;", "", ac)

    ## geneID
    eg <- item[grep("^DR\\s+GeneID", item)]
    eg <- gsub("\\D+", "", eg)

    ## embl
    emb <- item[grep("^DR\\s+EMBL", item)]
    emb <- sapply(emb, function(i) gsub("\\s+", "", strsplit(i, split=";")[[1]][3]))


    ## gene name
    gn <- item[grep("^GN", item)]
    gn <- gn[grep("Name", gn)]
    gn <- gsub("^GN\\s+Name=", "", gn)
    gn <- gsub("^GN\\s+OrderedLocusNames=", "", gn)
    gn <- gsub(";.*", "", gn)


    ## gene description
    de <- item[grep("DE\\s+RecName", item)][1]
    de <- gsub("^DE\\s[^=]+=", "", de)
    de <- gsub(";", "", de)

    ## GO
    go <- item[grep("^DR\\s+GO", item)]
    goid <- gsub("\\D+", "", go)
    if ( length(goid) != 0) {
        go <- paste("GO:", goid, sep="")
    } else {
        go <- NA
    }

    res <- data.frame(Uniprot=ac[1],
                      GeneID=eg[1],
                      EMBL=emb[1],
                      GeneName=gn[1],
                      FullName=de,
                      GO=go,
                      stringsAsFactors = FALSE,
                      check.names=FALSE)
    return(res)
}

