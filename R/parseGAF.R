##' parse GAF file
##'
##' given a GAF file, this function extracts the information from it 
##' @title parse.GAF
##' @param gafFile GAF file
##' @param nrows   a parameter
##' @return a list with two dataframes
##' @export
##' @import GO.db
##' @importFrom utils read.delim
##' @importFrom stats na.omit
##' @importFrom AnnotationDbi columns
parse.GAF <- function(gafFile,nrows=-1){
  gafFile <-read.GAF(gafFile)
  extr.gafFile <- gafFile[,c( "DB_Object_ID","GOID")]
  requireNamespace('GO.db')
  offspring.CC <- as.data.frame(GOCCOFFSPRING)
  offspring.MF <- as.data.frame(GOMFOFFSPRING)
  offspring.BP <- as.data.frame(GOBPOFFSPRING)
  offspring.total <- rbind(offspring.CC,offspring.MF,offspring.BP)
  names(offspring.total) <- c("OFFID","AncestorID")
  extr.offspring <- offspring.total[offspring.total$OFFID%in%extr.gafFile$GOID,]
  mer.info <- merge(extr.offspring,extr.gafFile,by.x="OFFID",by.y="GOID",all=T)
  mer.info <-na.omit(mer.info)
  mer.info <- mer.info[!duplicated(mer.info),]
  need.info1 <- mer.info[,c(2,3)]
  need.info2 <- mer.info[,c(1,3)]
  names(need.info1) <- c("GOID","DB_Object_ID")
  names(need.info2) <- c("GOID","DB_Object_ID")
  bind.info <- rbind(need.info1,need.info2)
  end.info <- bind.info[!duplicated(bind.info),]
  end.info <- end.info[order(end.info$GOID,end.info$DB_Object_ID),]
  go.ALL <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
  need.anno <- go.ALL[go.ALL$GOID%in%unique(end.info$GOID),]
  list(TERM2GENE = end.info[,c("GOID","DB_Object_ID")], TERM2NAME = need.anno[,c("GOID","TERM")])
  
}

##' @importFrom utils read.delim
read.GAF <- function(gafFile, nrows = -1) {
  cat("Reading ", gafFile, ": ", sep="")
  gafFile <- read.delim(gafFile, sep="\t", as.is=TRUE, quote="\"", fill=TRUE,
                        header=FALSE, nrows=nrows,comment.char="!")
  names(gafFile) = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GOID",
                     "DB:Reference(|DB:Reference)", "Evidence_Code", "With(or)Form", "Aspect","DB_Object_Name","DB_Object_Synonym(|Synonym)",
                     "DB_Object_Type","Taxon(|taxon)",
                     "Date","Assigned_by","Annotation_Extention","Gene_Product_Form_ID")
  
  cat("found", nrow(gafFile), "rows with classes:",
      paste(sapply(gafFile, class), collapse=", "), "\n")
  return(gafFile)
}

