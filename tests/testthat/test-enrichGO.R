library(clusterProfiler)

context("enrichGO")

test_that("non-ENTREZID universe is converted together with gene", {
    skip_if_not_installed("org.Hs.eg.db")

    raw_accnum_gene <- c(
        "O00411", "Q9BQP7", "Q8TCS8", "Q8WV74", "O14734", "Q96C36",
        "P48728", "Q5T440", "Q9NYK5", "Q92665", "Q92947", "Q969S9",
        "Q8WWH5", "Q9BSH4", "M0QWZ7", "P09543", "A0A0G2JL52", "Q7L8L6",
        "Q92506", "Q96GW9", "A0A0C4DGA2", "Q9BTZ2", "Q9NUB1", "P35914",
        "F6T1Q0", "P14174", "Q15067", "P30085", "J3KS15", "Q9BY49",
        "Q96GK7", "Q08426", "P52948", "A0A182DWF2", "O00161", "Q9P0Z9",
        "P22570", "P62760", "P43155", "Q99714", "P43897", "O75439",
        "Q96I99", "Q86X76", "F5H5I6", "Q8NE62", "H7BXY3", "Q9NVV4",
        "A0A0B4J1R2", "Q02252", "Q9BQ69", "P23378", "Q5JTZ9", "C9JQS9",
        "P04179", "P51649", "G3V325", "A0A2R8Y6Y7", "C9JRZ6", "P51659",
        "P36551", "Q9HCC0", "P04040", "P42765", "O15118", "F5GX62",
        "Q9NPJ3", "P55084", "O14561", "Q13011", "P35270", "P36957",
        "F6TLX2", "R4GN18", "Q9HBH1", "Q9Y624", "P53999", "Q9UIJ7",
        "Q8N5N7", "P24752", "Q6YN16", "P30084", "P23141", "P11413",
        "Q9BQA1", "H0YFD6", "P42126", "Q9UMS0", "Q13268", "Q00059",
        "Q86TX2", "P49748", "Q9HAV7", "O00330", "P61604", "P00390",
        "Q04837", "O76021", "P00352", "P27695", "Q9Y4W6", "P34897",
        "P28838", "P45954", "Q8TD30"
    )

    accnum_map <- suppressWarnings(bitr(
        raw_accnum_gene,
        fromType = "ACCNUM",
        toType = "ENTREZID",
        OrgDb = "org.Hs.eg.db"
    ))
    accnum_gene <- unique(accnum_map$ACCNUM)
    entrez_gene <- unique(accnum_map$ENTREZID)

    accnum_res <- enrichGO(
        accnum_gene,
        OrgDb = "org.Hs.eg.db",
        keyType = "ACCNUM",
        ont = "BP",
        universe = accnum_gene,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 1
    )
    expect_s4_class(accnum_res, "enrichResult")

    entrez_res <- enrichGO(
        entrez_gene,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENTREZID",
        ont = "BP",
        universe = entrez_gene,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 1
    )
    expect_s4_class(entrez_res, "enrichResult")

    compare_cols <- c(
        "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
        "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "Count"
    )

    expect_equal(
        as.data.frame(accnum_res)[, compare_cols],
        as.data.frame(entrez_res)[, compare_cols],
        tolerance = 1e-12
    )
})
