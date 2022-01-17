##' parse GAF files
##'
##' given a GAF file, this function extracts the information from it 
##' @title ParseGafFile
##' @param GafFile GAF file
##' @param nrows   a parameter
##' @return a list with two dataframes
##' @export
##' @importFrom GO.db GO.db
##' @importFrom GO.db GOCCOFFSPRING
##' @importFrom GO.db GOMFOFFSPRING
##' @importFrom GO.db GOBPOFFSPRING
##' @importFrom utils read.delim
##' @importFrom stats na.omit
##' @importFrom AnnotationDbi columns
ParseGafFile <- function(GafFile, nrows=-1) {
  
    GafFile <- ReadGafFile(GafFile)
    extr.gafFile <- GafFile[, c("DB_Object_ID", "GOID")]
  
    offspring.CC <- as.data.frame(GOCCOFFSPRING)
    offspring.MF <- as.data.frame(GOMFOFFSPRING)
    offspring.BP <- as.data.frame(GOBPOFFSPRING)
  
    offspring.total <- rbind(offspring.CC, offspring.MF, offspring.BP)
    names(offspring.total) <- c("OFFID", "AncestorID")
  
    extr.offspring <-
        offspring.total[offspring.total$OFFID %in% extr.gafFile$GOID, ]
  
    mer.info <-
        merge(
            extr.offspring,
            extr.gafFile,
            by.x = "OFFID",
            by.y = "GOID",
            all = T
            )
    mer.info <- na.omit(mer.info)
    mer.info <- mer.info[!duplicated(mer.info), ]
  
    need.info1 <- mer.info[, c(2, 3)]
    need.info2 <- mer.info[, c(1, 3)]
  
    names(need.info1) <- c("GOID", "DB_Object_ID")
    names(need.info2) <- c("GOID", "DB_Object_ID")
    bind.info <- rbind(need.info1, need.info2)
  
    end.info <- bind.info[!duplicated(bind.info),]
    end.info <- end.info[order(end.info$GOID, end.info$DB_Object_ID), ]
  
    go.ALL <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
    need.anno <- go.ALL[go.ALL$GOID %in% unique(end.info$GOID), ]
  
    list(TERM2GENE = end.info[, c("GOID", "DB_Object_ID")],  
         TERM2NAME = need.anno[, c("GOID", "TERM")])
  
}

##' @importFrom utils read.delim
ReadGafFile <- function(GafFile, nrows=-1) {
  
    cat("Reading ", GafFile, ": ", sep = "")
    GafFile <-
        read.delim(
            GafFile,
            sep = "\t",
            as.is = TRUE,
            quote = "\"",
            fill = TRUE,
            header = FALSE,
            nrows = nrows,
            comment.char = "!"
            )
    GafFile <- GafFile[, c(2, 3, 5, 7, 9, 10)]
    names(GafFile) <- c(
        "DB_Object_ID",
        "DB_Object_Symbol",
        "GOID",
        "Evidence_Code",
        "Aspect",
        "DB_Object_Name"
        )
    cat("found", nrow(GafFile), "rows in this GAF file")
    return(GafFile)
}

