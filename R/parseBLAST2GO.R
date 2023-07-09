##' parse BLAST2GO files
##'
##' given a BLAST2GO file, this function extracts the information from it and make it use for TERM2GENE.
##' @title parse_blast2go
##' @param Blast2GOFile BLAST2GO file
##' @return a data frame with two columns: geneID and GO.ID
##' @export
##' @importFrom GO.db GO.db
##' @importFrom utils read.delim
##' @importFrom AnnotationDbi columns
parse_blast2go <- function(Blast2GOFile) {
    blast2go <- read.table(Blast2GOFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    go_annotation_data <- blast2go[, c("Sequence.Name", "Blast.Top.Hit.GOs")]
    go_annotation_data <- tidyr::separate_rows(go_annotation_data, `Blast.Top.Hit.GOs`, sep = ", ")
    names(go_annotation_data) <- c("geneID", "GO.ID")
    return(go_annotation_data)
}

