##' new parse GAF files
##'
##' given a GAF file, this function extracts the information from it 
##' @title parse_gff
##' @param GafFile GAF file
##' @param nrows   a parameter
##' @return a list with two dataframes
##' @export
##' @importFrom utils read.delim
##' @importFrom GO.db GO.db
##' @importFrom AnnotationDbi columns



parse_gff <- function(GafFile, nrows = -1) {
  GafFile <- read.gff(GafFile)
  new.data.frame <- GafFile[, c("GOID", "DB_Object_ID")]
  
  ##use buildGOmap function to get information related
  build.df <- buildGOmap(new.data.frame)
  
  ##rename this colnames to be same with the buildGOmap function to facilatate bind information.
  names(new.data.frame) <- c("GO", "Gene")
  
  ##bind the information needed
  bind.info <- rbind(build.df, new.data.frame)
  
  ##process  data
  bind.info <- bind.info[order(bind.info$GO, bind.info$Gene),]
  bind.info <- bind.info[!duplicated(bind.info),]
  bind.info[, "Gene"] <- as.character(bind.info$Gene)
  
  go.ALL <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
  need.anno <- go.ALL[go.ALL$GOID %in% unique(bind.info$GO),]
  
  list(TERM2GENE = bind.info[, c("GO", "Gene")],
       TERM2NAME = need.anno[, c("GOID", "TERM")])
  
}
read.gff <- function(GafFile, nrows=-1) {
  
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
setwd("D:/")
