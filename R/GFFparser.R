#' read GFF file and build gene information table
#'
#' given a GFF file, this function extracts information from it and save it in working directory
#' @title Gff2GeneTable
#' @param gffFile GFF file
#' @param compress compress file or not
#' @return file save.
#' @export
#' @author Yu Guangchuang
Gff2GeneTable <- function(gffFile, compress=TRUE) {
    gff <- readGff(gffFile)

    ## NCBI/RefSeq GFF3 writes `GeneID=...` while Ensembl (and GTF-derived GFF3)
    ## writes `gene_id=...`. Looking only for `GeneID` on an Ensembl file leaves
    ## every key NA, and merging two all-NA key columns is a cartesian product:
    ## silently useless on a small file, and on a real genome it overflows with
    ## "negative length vectors are not allowed" (#193).
    gid_field <- detect_gff_field(gff$attributes, c("GeneID", "gene_id"))
    if (is.null(gid_field)) {
        stop("no gene identifier found in the GFF attributes; ",
             "Gff2GeneTable() expects a `GeneID=` attribute (NCBI/RefSeq style) ",
             "or a `gene_id=` attribute (Ensembl/GTF style)")
    }

    GeneID <- data.frame(GeneID=getGffAttribution(gff$attributes, field=gid_field)
                         )
    ## GI2GeneID <- data.frame(GI=getGffAttribution(gff$attributes, field="GI"),
    ##                        GeneID=getGffAttribution(gff$attributes, field="GeneID")
    ##                                    #,
    ##                                    #Product=getGffAttribution(gff$attributes, field="product")
    ##                        )
    ## GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$GI),]
    ## GI2GeneID <- GI2GeneID[!is.na(GI2GeneID$Gene),]

    geneInfo <- gff[gff$feature == "gene",]
    if (nrow(geneInfo) == 0) {
        stop("the GFF file has no rows with feature == \"gene\"")
    }
    geneInfo <- geneInfo[, c("seqname", "start", "end", "strand", "attributes")]
    geneInfo$GeneID <- getGffAttribution(geneInfo$attributes, field=gid_field)

    ## the gene symbol is `gene=` in NCBI files and `Name=` in Ensembl files
    gname_field <- detect_gff_field(geneInfo$attributes, c("gene", "Name", "gene_name"))
    geneInfo$GeneName <- if (is.null(gname_field)) {
        NA_character_
    } else {
        getGffAttribution(geneInfo$attributes, field=gname_field)
    }
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

#' Return the first GFF attribute name that is actually populated
#'
#' Tries each candidate in turn and returns the first that appears as an exact
#' attribute key in the file *and* yields at least one non-empty value, or NULL
#' when none of them do.
#'
#' The exact-key test matters: `getGffAttribution()` has a heuristic fallback
#' that also looks inside values, and on an Ensembl record it happily returns
#' the identifier out of `ID=gene:ENSG...` when asked for `gene`. Requiring the
#' key to exist keeps `gene=` (NCBI symbol) from shadowing `Name=` (Ensembl
#' symbol).
#'
#' @param attributes the attributes column of a parsed GFF
#' @param candidates attribute names to try, in order of preference
#' @return a single string, or NULL
#' @noRd
detect_gff_field <- function(attributes, candidates) {
    for (field in candidates) {
        if (!any(grepl(paste0("(^|;)", field, "="), attributes))) {
            next
        }
        value <- getGffAttribution(attributes, field = field)
        if (any(!is.na(value) & nzchar(value))) {
            return(field)
        }
    }
    NULL
}

##
## M5005 GFF file was downloaded from:
## ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Streptococcus_pyogenes_MGAS5005_uid58337/
##
##
## Gff2GeneTable("NC_007297.gff")
##
##
#' @importFrom utils read.table
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
