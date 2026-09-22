## Gff2GeneTable(): NCBI/RefSeq and Ensembl attribute styles (#193)

write_gff <- function(lines, dir) {
    path <- file.path(dir, "test.gff3")
    writeLines(lines, path)
    path
}

## Gff2GeneTable() saves geneTable.rda into the working directory, so run it in a
## scratch dir and read the table back from there.
run_gff <- function(gff, dir) {
    old <- setwd(dir)
    on.exit(setwd(old), add = TRUE)
    Gff2GeneTable(gff, compress = FALSE)
    e <- new.env()
    load(file.path(dir, "geneTable.rda"), envir = e)
    get("geneTable", envir = e)
}

ensembl_lines <- c(
    "##gff-version 3",
    "1\tensembl\tgene\t1000\t1100\t.\t+\t.\tID=gene:ENSFAKE00000000001;Name=GENEA;gene_id=ENSFAKE00000000001;gene_biotype=protein_coding",
    "1\tensembl\tCDS\t1000\t1000\t.\t+\t0\tID=CDS:ENSFAKE00000000001.1;Parent=transcript:ENSFAKE00000000001.1;gene_id=ENSFAKE00000000001",
    "1\tensembl\tgene\t2000\t2100\t.\t-\t.\tID=gene:ENSFAKE00000000002;Name=GENEB;gene_id=ENSFAKE00000000002;gene_biotype=protein_coding",
    "1\tensembl\tCDS\t2000\t2000\t.\t-\t0\tID=CDS:ENSFAKE00000000002.1;Parent=transcript:ENSFAKE00000000002.1;gene_id=ENSFAKE00000000002"
)

ncbi_lines <- c(
    "##gff-version 3",
    "NC_1\tRefSeq\tgene\t1000\t1100\t.\t+\t.\tGeneID=100001;gene=SYMA",
    "NC_1\tRefSeq\tCDS\t1000\t1000\t.\t+\t0\tGeneID=100001;product=hypothetical",
    "NC_1\tRefSeq\tgene\t2000\t2100\t.\t-\t.\tGeneID=100002;gene=SYMB",
    "NC_1\tRefSeq\tCDS\t2000\t2000\t.\t-\t0\tGeneID=100002;product=hypothetical"
)

test_that("Gff2GeneTable() reads Ensembl-style gene_id attributes (#193)", {
    dir <- tempfile("gffens")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    gff <- write_gff(ensembl_lines, dir)

    gt <- run_gff(gff, dir)

    expect_equal(nrow(gt), 2)
    expect_false(any(is.na(gt$GeneID)))
    expect_setequal(gt$GeneID,
        c("ENSFAKE00000000001", "ENSFAKE00000000002"))
    # the symbol comes from Name=, not from the gene:... part of the ID
    expect_setequal(gt$GeneName, c("GENEA", "GENEB"))
    expect_true(all(gt$start %in% c(1000, 2000)))
})

test_that("Gff2GeneTable() still reads NCBI-style GeneID/gene attributes", {
    dir <- tempfile("gffncbi")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    gff <- write_gff(ncbi_lines, dir)

    gt <- run_gff(gff, dir)

    expect_equal(nrow(gt), 2)
    expect_setequal(gt$GeneID, c("100001", "100002"))
    expect_setequal(gt$GeneName, c("SYMA", "SYMB"))
})

test_that("Gff2GeneTable() reports an unrecognised attribute style", {
    dir <- tempfile("gffweird")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    gff <- write_gff(c("##gff-version 3",
                       "1\tx\tgene\t1\t9\t.\t+\t.\tfoo=bar;baz=qux"), dir)

    old <- setwd(dir)
    on.exit(setwd(old), add = TRUE)
    expect_error(Gff2GeneTable(gff, compress = FALSE), "no gene identifier")
})

test_that("detect_gff_field() requires an exact attribute key", {
    ens <- "ID=gene:ENSG1;Name=SYM;gene_id=ENSG1"
    expect_equal(clusterProfiler:::detect_gff_field(ens, c("GeneID", "gene_id")), "gene_id")
    # `gene` must not match the `gene:` inside ID=gene:...
    expect_equal(clusterProfiler:::detect_gff_field(ens, c("gene", "Name", "gene_name")), "Name")

    ncbi <- "GeneID=1;gene=SYM"
    expect_equal(clusterProfiler:::detect_gff_field(ncbi, c("GeneID", "gene_id")), "GeneID")
    expect_equal(clusterProfiler:::detect_gff_field(ncbi, c("gene", "Name")), "gene")

    expect_null(clusterProfiler:::detect_gff_field("foo=bar", c("GeneID", "gene_id")))
})
