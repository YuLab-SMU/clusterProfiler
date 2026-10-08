## Pathway Commons gene sets whose name field holds the "; " separator (#789)

## A Pathway Commons GMT line is the identifier of the gene set, a description
## made of labelled fields, and then the genes. The identifier is a plain
## compact URI, or an http address whose last path segment is the identifier.
pc_line <- function(id, description, genes) {
    paste(c(id, description, genes), collapse = "\t")
}

pc_gmt <- function(lines) {
    path <- tempfile(fileext = ".gmt")
    writeLines(lines, path)
    path
}

test_that("read.gmt.pc() keeps the fields aligned when a name contains '; ' (#789)", {
    gmt <- pc_gmt(c(
        pc_line("reactome:R-HSA-5654221",
                "name: Phospholipase C-mediated cascade; FGFR2; datasource: reactome; organism: 9606; idtype: hgnc.symbol",
                c("FGFR2", "PLCG1")),
        pc_line("http://bioregistry.io/kegg.pathway:hsa00010",
                "name: Glycolysis / Gluconeogenesis; datasource: kegg; organism: 9606; idtype: hgnc.symbol",
                "ADORA1"),
        pc_line("biofactoid:00e3ac12-8c96-45df-a59c-e212724962b6",
                "name: Arginine: Glycine Amidinotransferase Deficiency (AGAT Deficiency); datasource: pathbank; organism: 9606; idtype: hgnc.symbol",
                "AGAT")
    ))

    x <- expect_silent(read.gmt.pc(gmt))

    ## an address is shortened to the identifier it ends with
    expect_equal(unique(x$id),
                 c("reactome:R-HSA-5654221", "kegg.pathway:hsa00010",
                   "biofactoid:00e3ac12-8c96-45df-a59c-e212724962b6"))
    ## a '; ' inside the name is part of the name, a colon is no separator either
    expect_equal(x$name[1], "Phospholipase C-mediated cascade; FGFR2")
    expect_equal(x$name[4],
                 "Arginine: Glycine Amidinotransferase Deficiency (AGAT Deficiency)")
    expect_equal(unique(x$datasource), c("reactome", "kegg", "pathbank"))
    expect_equal(unique(x$organism), "9606")
    expect_equal(unique(x$idtype), "hgnc.symbol")
})

test_that("read.gmt.pc() reads a description that declares fewer fields (#789)", {
    gmt <- pc_gmt(c(
        pc_line("wp:WP100", "name: Cell cycle; organism: 9606; idtype: hgnc.symbol",
                c("CDK1", "CCNA2")),
        pc_line("wp:WP106",
                "name: Programmed cell death; datasource: wikipathways; organism: 9606; idtype: hgnc.symbol",
                "CASP8")
    ))

    x <- expect_silent(read.gmt.pc(gmt))

    expect_equal(x$datasource, c(NA, NA, "wikipathways"))
    expect_equal(unique(x$organism), "9606")
    expect_equal(x$name, c("Cell cycle", "Cell cycle", "Programmed cell death"))

    y <- expect_silent(read.gmt.pc(gmt, output = "GSON"))
    expect_s4_class(y, "GSON")
    expect_equal(y@species, "Homo sapiens")
})

test_that("read.gmt.pc() takes the species from a gene set that declares one (#789)", {
    gmt <- pc_gmt(c(
        pc_line("netpath5_pathways:NP_001625",
                "name: DNA double-strand break signaling; idtype: hgnc.symbol", "ATM"),
        pc_line("reactome:R-HSA-5654227",
                "name: Phospholipase C-mediated cascade; FGFR3; datasource: reactome; organism: 9606; idtype: hgnc.symbol",
                c("FGFR3", "PLCG1"))
    ))

    y <- expect_silent(read.gmt.pc(gmt, output = "gson"))
    expect_equal(y@species, "Homo sapiens")
    expect_equal(y@gsid2name$name[y@gsid2name$gsid == "reactome:R-HSA-5654227"],
                 "Phospholipase C-mediated cascade; FGFR3")
})
