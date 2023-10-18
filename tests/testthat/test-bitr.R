library(clusterProfiler)

context("bitr")

test_that("bitr", {
    expect_true('P22223' %in% bitr('1001', 'ENTREZID', 'UNIPROT', 'org.Hs.eg.db')[,2])
    ## maybe network problem
    res <- tryCatch(bitr_kegg('1001', 'ncbi-geneid', 'uniprot', 'hsa'), error = function(e) NULL)
    if (!is.null(res)) expect_true('P22223' %in% res[,2])
})

