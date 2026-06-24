library(clusterProfiler)

context("bitr")

test_that("bitr", {
    expect_true('P22223' %in% bitr('1001', 'ENTREZID', 'UNIPROT', 'org.Hs.eg.db')[,2])
    ## maybe network problem
    res <- suppressWarnings(
        tryCatch(bitr_kegg('1001', 'ncbi-geneid', 'uniprot', 'hsa'), error = function(e) NULL)
    )
    if (!is.null(res)) expect_true('P22223' %in% res[,2])
})

test_that("bitr_kegg supports UniProt to KO conversion", {
    local_mocked_bindings(
        kegg_rest = function(rest_url) {
            if (grepl("/conv/uniprot/hsa$", rest_url)) {
                return(data.frame(
                    from = c("hsa:10458", "hsa:3630"),
                    to = c("up:P31946", "up:P01308"),
                    stringsAsFactors = FALSE
                ))
            }
            if (grepl("/link/ko/hsa$", rest_url)) {
                return(data.frame(
                    from = c("hsa:10458", "hsa:3630"),
                    to = c("ko:K04456", "ko:K04526"),
                    stringsAsFactors = FALSE
                ))
            }
            stop("unexpected url: ", rest_url)
        },
        .package = "clusterProfiler"
    )

    res <- bitr_kegg(c("P31946", "P01308"), "uniprot", "ko", "hsa")

    expect_equal(res$uniprot, c("P31946", "P01308"))
    expect_equal(res$ko, c("K04456", "K04526"))
})

test_that("bitr_kegg supports KO to external ID conversion", {
    local_mocked_bindings(
        kegg_rest = function(rest_url) {
            if (grepl("/link/hsa/ko$", rest_url)) {
                return(data.frame(
                    from = c("ko:K04456", "ko:K04526"),
                    to = c("hsa:10458", "hsa:3630"),
                    stringsAsFactors = FALSE
                ))
            }
            if (grepl("/conv/uniprot/hsa$", rest_url)) {
                return(data.frame(
                    from = c("hsa:10458", "hsa:3630"),
                    to = c("up:P31946", "up:P01308"),
                    stringsAsFactors = FALSE
                ))
            }
            if (grepl("/conv/ncbi-geneid/hsa$", rest_url)) {
                return(data.frame(
                    from = c("hsa:10458", "hsa:3630"),
                    to = c("ncbi-geneid:10458", "ncbi-geneid:3630"),
                    stringsAsFactors = FALSE
                ))
            }
            stop("unexpected url: ", rest_url)
        },
        .package = "clusterProfiler"
    )

    ko2uniprot <- bitr_kegg(c("K04456", "K04526"), "ko", "uniprot", "hsa")
    ko2geneid <- bitr_kegg(c("K04456", "K04526"), "ko", "ncbi-geneid", "hsa")

    expect_equal(ko2uniprot$ko, c("K04456", "K04526"))
    expect_equal(ko2uniprot$uniprot, c("P31946", "P01308"))
    expect_equal(ko2geneid$ko, c("K04456", "K04526"))
    expect_equal(ko2geneid$`ncbi-geneid`, c("10458", "3630"))
})

test_that("bitr_kegg supports KO as KEGG keyType for pathway conversion", {
    local_mocked_bindings(
        download_KEGG = function(species, keggType = "KEGG", keyType = "kegg") {
            expect_equal(species, "hsa")
            expect_equal(keggType, "KEGG")
            expect_equal(keyType, "ko")
            list(
                KEGGPATHID2EXTID = data.frame(
                    from = c("hsa04010", "hsa04010", "hsa04910"),
                    to = c("K04456", "K04526", "K04526"),
                    stringsAsFactors = FALSE
                ),
                KEGGPATHID2NAME = data.frame(
                    from = c("hsa04010", "hsa04910"),
                    to = c("MAPK signaling pathway - Homo sapiens", "Insulin signaling pathway - Homo sapiens"),
                    stringsAsFactors = FALSE
                )
            )
        },
        .package = "clusterProfiler"
    )

    res <- bitr_kegg("K04526", "ko", "Path", "hsa")

    expect_equal(res$ko, c("K04526", "K04526"))
    expect_equal(res$Path, c("hsa04010", "hsa04910"))
})

test_that("gson_KO is exported and builds a KO GSON object", {
    local_mocked_bindings(
        kegg_rest = function(rest_url) {
            if (grepl("/link/ko/pathway$", rest_url)) {
                return(data.frame(
                    from = c("path:map00010", "path:map00020"),
                    to = c("ko:K00844", "ko:K00116"),
                    stringsAsFactors = FALSE
                ))
            }
            if (grepl("/list/pathway$", rest_url)) {
                return(data.frame(
                    from = c("path:map00010", "path:map00020"),
                    to = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)"),
                    stringsAsFactors = FALSE
                ))
            }
            stop("unexpected url: ", rest_url)
        },
        .package = "clusterProfiler"
    )
    local_mocked_bindings(
        yread = function(url) c("ko             Release 116.0+/06-24, Jun 26"),
        .package = "yulab.utils"
    )

    x <- gson_KO()

    expect_s4_class(x, "GSON")
    expect_equal(x@species, "KEGG Orthology")
    expect_equal(x@keytype, "kegg_orthology")
    expect_true(all(c("K00844", "K00116") %in% x@gsid2gene$gene))
})

