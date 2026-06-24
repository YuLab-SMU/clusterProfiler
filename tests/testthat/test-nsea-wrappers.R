library(clusterProfiler)

context("nsea wrappers")

toy_kegg_gson <- function() {
    gson::gson(
        gsid2gene = data.frame(
            gsid = c("hsa00010", "hsa00010", "hsa00020", "hsa00020"),
            gene = c("g1", "g2", "g3", "g4"),
            stringsAsFactors = FALSE
        ),
        gsid2name = data.frame(
            gsid = c("hsa00010", "hsa00020"),
            name = c("Toy glycolysis", "Toy citrate cycle"),
            stringsAsFactors = FALSE
        ),
        species = "hsa",
        gsname = "KEGG",
        keytype = "kegg",
        version = "toy",
        accessed_date = as.character(Sys.Date())
    )
}

toy_mkegg_gson <- function() {
    gson::gson(
        gsid2gene = data.frame(
            gsid = c("M00001", "M00001", "M00002", "M00002"),
            gene = c("g1", "g2", "g3", "g4"),
            stringsAsFactors = FALSE
        ),
        gsid2name = data.frame(
            gsid = c("M00001", "M00002"),
            name = c("Toy module 1", "Toy module 2"),
            stringsAsFactors = FALSE
        ),
        species = "hsa",
        gsname = "MKEGG",
        keytype = "kegg",
        version = "toy",
        accessed_date = as.character(Sys.Date())
    )
}

toy_wp_gson <- function() {
    gson::gson(
        gsid2gene = data.frame(
            gsid = c("WP1", "WP1", "WP2", "WP2"),
            gene = c("g1", "g2", "g3", "g4"),
            stringsAsFactors = FALSE
        ),
        gsid2name = data.frame(
            gsid = c("WP1", "WP2"),
            name = c("Toy pathway 1", "Toy pathway 2"),
            stringsAsFactors = FALSE
        ),
        species = "Homo sapiens",
        gsname = "WikiPathways",
        keytype = "ENTREZID",
        version = "toy",
        accessed_date = as.character(Sys.Date())
    )
}

toy_network <- data.frame(
    from = c("g1", "g2", "g3", "g4"),
    to = c("g2", "g3", "g4", "g1"),
    weight = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
)

toy_gene_list <- c(g1 = 4, g2 = 3, g3 = 2, g4 = 1)

toy_layer_networks <- list(
    transcriptome = toy_network,
    proteome = toy_network
)

toy_couplings <- data.frame(
    from_layer = c("transcriptome", "transcriptome", "proteome", "proteome"),
    from_id = c("g1", "g2", "g1", "g2"),
    to_layer = c("proteome", "proteome", "transcriptome", "transcriptome"),
    to_id = c("g1", "g2", "g1", "g2"),
    weight = 1,
    stringsAsFactors = FALSE
)

toy_seed_list <- list(
    transcriptome = toy_gene_list,
    proteome = c(g1 = 2, g2 = 1.5, g3 = 1, g4 = 0.5)
)

test_that("nseGO returns nseaResult and marks SYMBOL input readable", {
    skip_if_not_installed("org.Hs.eg.db")

    geneList <- c(
        TP53 = 5,
        EGFR = 4,
        VEGFA = 3,
        MYC = 2,
        AKT1 = 1.5,
        BRCA1 = 1
    )
    network <- data.frame(
        from = c("TP53", "EGFR", "VEGFA", "MYC", "AKT1", "BRCA1"),
        to = c("EGFR", "VEGFA", "MYC", "AKT1", "BRCA1", "TP53"),
        weight = 1,
        stringsAsFactors = FALSE
    )

    res <- nseGO(
        geneList = geneList,
        network = network,
        ont = "BP",
        OrgDb = "org.Hs.eg.db",
        keyType = "SYMBOL",
        minGSSize = 1,
        maxGSSize = 500,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "nseaResult")
    expect_true(res@readable)
    expect_equal(res@setType, "BP")
    expect_equal(res@keytype, "SYMBOL")
    expect_true(nrow(as.data.frame(res)) >= 1)
})

test_that("mnseGO returns mnseaResult with GO metadata", {
    skip_if_not_installed("org.Hs.eg.db")

    seed_list <- list(
        transcriptome = c(TP53 = 5, EGFR = 4, VEGFA = 3, MYC = 2),
        proteome = c(TP53 = 2.5, EGFR = 2, VEGFA = 1.5, MYC = 1)
    )
    networks <- list(
        transcriptome = data.frame(
            from = c("TP53", "EGFR", "VEGFA", "MYC"),
            to = c("EGFR", "VEGFA", "MYC", "TP53"),
            weight = 1,
            stringsAsFactors = FALSE
        ),
        proteome = data.frame(
            from = c("TP53", "EGFR", "VEGFA", "MYC"),
            to = c("EGFR", "VEGFA", "MYC", "TP53"),
            weight = 1,
            stringsAsFactors = FALSE
        )
    )
    couplings <- data.frame(
        from_layer = c("transcriptome", "transcriptome", "proteome", "proteome"),
        from_id = c("TP53", "EGFR", "TP53", "EGFR"),
        to_layer = c("proteome", "proteome", "transcriptome", "transcriptome"),
        to_id = c("TP53", "EGFR", "TP53", "EGFR"),
        weight = 1,
        stringsAsFactors = FALSE
    )

    res <- mnseGO(
        seed_list = seed_list,
        networks = networks,
        couplings = couplings,
        ont = "BP",
        OrgDb = "org.Hs.eg.db",
        keyType = "SYMBOL",
        minGSSize = 1,
        maxGSSize = 500,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "mnseaResult")
    expect_equal(res@setType, "BP")
    expect_equal(res@keytype, "SYMBOL")
    expect_true(res@readable)
    expect_true(length(res@layer_scores) == 2)
})

test_that("nseKEGG accepts GSON input and appends KEGG metadata", {
    local_mocked_bindings(append_kegg_category = function(x) x, .package = "clusterProfiler")

    res <- nseKEGG(
        geneList = toy_gene_list,
        network = toy_network,
        organism = toy_kegg_gson(),
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "nseaResult")
    expect_equal(res@setType, "KEGG")
    expect_equal(res@organism, "hsa")
    expect_equal(res@keytype, "kegg")
})

test_that("mnseKEGG accepts GSON input and returns mnseaResult", {
    local_mocked_bindings(append_kegg_category = function(x) x, .package = "clusterProfiler")

    res <- mnseKEGG(
        seed_list = toy_seed_list,
        networks = toy_layer_networks,
        couplings = toy_couplings,
        organism = toy_kegg_gson(),
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "mnseaResult")
    expect_equal(res@setType, "KEGG")
    expect_equal(res@organism, "hsa")
    expect_equal(res@keytype, "kegg")
    expect_true(length(res@layer_scores) == 2)
})

test_that("nseMKEGG accepts GSON input and appends module metadata", {
    res <- nseMKEGG(
        geneList = toy_gene_list,
        network = toy_network,
        organism = toy_mkegg_gson(),
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "nseaResult")
    expect_equal(res@setType, "MKEGG")
    expect_equal(res@organism, "hsa")
    expect_equal(res@keytype, "kegg")
})

test_that("mnseMKEGG accepts GSON input and returns mnseaResult", {
    res <- mnseMKEGG(
        seed_list = toy_seed_list,
        networks = toy_layer_networks,
        couplings = toy_couplings,
        organism = toy_mkegg_gson(),
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "mnseaResult")
    expect_equal(res@setType, "MKEGG")
    expect_equal(res@organism, "hsa")
    expect_equal(res@keytype, "kegg")
    expect_true(length(res@layer_scores) == 2)
})

test_that("nseWP accepts GSON input and appends WikiPathways metadata", {
    res <- nseWP(
        geneList = toy_gene_list,
        network = toy_network,
        organism = toy_wp_gson(),
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "nseaResult")
    expect_equal(res@setType, "WikiPathways")
    expect_equal(res@organism, "Homo sapiens")
    expect_equal(res@keytype, "ENTREZID")
})

test_that("mnseWP accepts GSON input and returns mnseaResult", {
    res <- mnseWP(
        seed_list = toy_seed_list,
        networks = toy_layer_networks,
        couplings = toy_couplings,
        organism = toy_wp_gson(),
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "permute",
        nPerm = 10,
        verbose = FALSE
    )

    expect_s4_class(res, "mnseaResult")
    expect_equal(res@setType, "WikiPathways")
    expect_equal(res@organism, "Homo sapiens")
    expect_equal(res@keytype, "ENTREZID")
    expect_true(length(res@layer_scores) == 2)
})
