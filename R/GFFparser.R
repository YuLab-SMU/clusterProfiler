##' read GFF file and build gene information table
##'
##' given a GFF file, this function extracts information from it and save it in working directory
##' @title Gff2GeneTable
##' @param gffFile GFF file
##' @param compress compress file or not
##' @return file save.
##' @export
##' @author Yu Guangchuang
Gff2GeneTable <- function(gffFile, compress=TRUE) {
    gff <- readGff(gffFile)

    GeneID <- data.frame(GeneID=getGffAttribution(gff$attributes, field="GeneID")
                         )
    ## GI2GeneID <- data.frame(GI=getGffAttribution(gff$attributes, field="GI"),
    ##                        GeneID=getGffAttribution(gff$attributes, field="GeneID")
    ##                                    #,
    ##                                    #Product=getGffAttribution(gff$attributes, field="product")
    ##                        )
    ## GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$GI),]
    ## GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$Gene),]

    geneInfo <- gff[gff$feature == "gene",]
    geneInfo <- geneInfo[, c("seqname", "start", "end", "strand", "attributes")]
    geneInfo$GeneID <- getGffAttribution(geneInfo$attributes, field="GeneID")
    geneInfo$GeneName <- getGffAttribution(geneInfo$attributes, field="gene")
    geneInfo$Locus <- getGffAttribution(geneInfo$attributes, field="locus_tag")
    geneInfo$GeneName[is.na(geneInfo$GeneName)] <- "-"

    geneInfo <- geneInfo[, -5] ## abondom "attributes" column.
    ## geneTable <- merge(GI2GeneID, geneInfo, by.x="GeneID", by.y="GeneID")
    geneTable <- merge(GeneID, geneInfo, by.x="GeneID", by.y="GeneID")
    geneTable <- unique(geneTable)
    if (compress) {
        save(geneTable, file="geneTable.rda", compress="xz")
    } else {
        save(geneTable, file="geneTable.rda")
    }
    print("Gene Table file save in the working directory.")
}

##
## M5005 GFF file was downloaded from:
## ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Streptococcus_pyogenes_MGAS5005_uid58337/
##
##
## Gff2GeneTable("NC_007297.gff")
##
##
##' @importFrom utils read.table
readGff <- function(gffFile, nrows = -1) {
    cat("Reading ", gffFile, ": ", sep="")
    gff <- read.table(gffFile, sep="\t", as.is=TRUE, quote="\"", fill=TRUE,
                      header=FALSE, comment.char="#", nrows=nrows,
                      colClasses=c("character", "character", "character", "integer",
                      "integer", "character", "character", "character", "character"))
    colnames(gff) = c("seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes")
    cat("found", nrow(gff), "rows with classes:",
        paste(sapply(gff, class), collapse=", "), "\n")
    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
    return(gff)
}

getGffAttribution <- function (x, field, attrsep = ";") {
    s = strsplit(x, split = attrsep, fixed = TRUE)
    sapply(s, function(atts) {
        a = strsplit(atts, split = "=", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
        } else {
            b = sapply(a, function(atts) {
                strsplit(atts[2], split = ",", fixed = TRUE)
            })

            rv = as.character(NA)
            sapply(b, function(atts) {
                secA <- strsplit(atts, split = ":", fixed = TRUE)
                m = match(field, sapply(secA, "[", 1))

                if (!is.na(m)) {
                    rv <<- secA[[m]][2]
                }

            })


        }
        return(rv)
    })
}
