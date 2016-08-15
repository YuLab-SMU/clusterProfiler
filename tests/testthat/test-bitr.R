library(clusterProfiler)

context("bitr")

test_that("bitr", {
    expect_true('P22223' %in% bitr('1001', 'ENTREZID', 'UNIPROT', 'org.Hs.eg.db')[,2])
    expect_true('P22223' %in% bitr_kegg('1001', 'ncbi-geneid', 'uniprot', 'hsa')[,2])
})

