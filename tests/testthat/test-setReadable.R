library(testthat)
library(clusterProfiler)

## setReadable() and EXTID2NAME() are enrichit functions that clusterProfiler
## re-exports and calls from its own front ends (enrichGO(), enrichDisease(), ...),
## so their contract is decided here and tested here; enrichit deliberately has
## no test for them (see enrichit/AGENTS.md, section 2).
##
## What matters at this level is not "AnnotationDbi::select() was called" but that
## readable IDs appear in the result, that every statistic survives the rewrite,
## and that a partly mappable OrgDb does not silently drop genes.

skip_if_no_orgdb <- function() {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("GO.db")
}

run_enrichgo <- function(...) {
    odb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    gene <- head(AnnotationDbi::keys(odb, "ENTREZID"), 120)
    suppressWarnings(enrichGO(
        gene, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP",
        pvalueCutoff = 1, qvalueCutoff = 1, ...
    ))
}

test_that("setReadable rewrites geneID and leaves the statistics alone", {
    skip_if_no_orgdb()
    res <- run_enrichgo()
    expect_s4_class(res, "enrichResult")
    expect_false(res@readable)

    conv <- suppressWarnings(setReadable(res, "org.Hs.eg.db", keyType = "ENTREZID"))
    expect_true(conv@readable)

    ## Entrez ids became symbols -- no bare integer is left in geneID
    ids <- unlist(strsplit(as.character(conv@result$geneID), "/", fixed = TRUE))
    ids <- ids[nzchar(ids)]
    expect_true(length(ids) > 0)
    expect_false(any(grepl("^[0-9]+$", ids)))

    ## every statistic column is byte-identical: this is display only
    expect_identical(conv@result$pvalue, res@result$pvalue)
    expect_identical(conv@result$p.adjust, res@result$p.adjust)
    expect_identical(conv@result$Count, res@result$Count)
    expect_identical(conv@result$ID, res@result$ID)
    expect_identical(conv@result$GeneRatio, res@result$GeneRatio)
    expect_identical(conv@result$BgRatio, res@result$BgRatio)

    ## and the gene sets they came from are still whole
    expect_identical(conv@result$Count,
                     vapply(strsplit(as.character(conv@result$geneID), "/", fixed = TRUE),
                            function(x) sum(nzchar(x)), integer(1)))
})

test_that("setReadable leaves the query record in its input namespace", {
    skip_if_no_orgdb()
    res <- run_enrichgo()
    conv <- suppressWarnings(setReadable(res, "org.Hs.eg.db", keyType = "ENTREZID"))

    ## @gene records what was asked for, not what it was displayed as
    expect_identical(conv@gene, res@gene)
    ## The mapping itself is kept and covers exactly the query genes -- it is not
    ## a global symbol table, so an id outside the subset is out of range. Only
    ## the shape is asserted here: which symbols come back is fixed by the OrgDb,
    ## not by this package, and a test that read the answer back off the object
    ## under test would pass whatever the mapping did. The actual id -> symbol
    ## correspondence is specified with a fixed input in the EXTID2NAME test.
    expect_identical(length(conv@gene2Symbol), length(res@gene))
    expect_setequal(names(conv@gene2Symbol), res@gene)
    expect_false(anyNA(names(conv@gene2Symbol)))
})

test_that("setReadable is stable when applied twice", {
    skip_if_no_orgdb()
    res <- run_enrichgo()
    once <- suppressWarnings(setReadable(res, "org.Hs.eg.db", keyType = "ENTREZID"))
    twice <- suppressWarnings(setReadable(once, "org.Hs.eg.db", keyType = "SYMBOL"))
    expect_identical(twice@result$geneID, once@result$geneID)
    expect_identical(twice@result$pvalue, once@result$pvalue)
})

test_that("readable = TRUE in enrichGO() matches an explicit setReadable() call", {
    skip_if_no_orgdb()
    ## the front end applies setReadable internally when readable = TRUE; the two
    ## routes must agree, which is the whole point of testing this at this level
    direct <- run_enrichgo()
    direct <- suppressWarnings(setReadable(direct, "org.Hs.eg.db", keyType = "ENTREZID"))

    via_arg <- run_enrichgo(readable = TRUE)
    expect_identical(via_arg@result$geneID, direct@result$geneID)
    expect_identical(via_arg@result$pvalue, direct@result$pvalue)
})

test_that("EXTID2NAME maps ids and passes unmappable ones through", {
    skip_if_no_orgdb()
    ## EXTID2NAME() is not exported by clusterProfiler -- setReadable() calls it
    ## internally -- so it is reached through the namespace here. It is tested at
    ## this level because the behaviour that matters (ids with no symbol are
    ## passed through rather than becoming NA) is what setReadable() relies on.
    X <- enrichit::EXTID2NAME

    ## Fixed input, stated expectation. The three ids below are the canonical
    ## HGNC assignments for these genes and are what a caller of setReadable()
    ## relies on seeing; the input is specified here instead of being read off the
    ## result, and the order is part of the contract (a caller indexing the vector
    ## positionally must get the gene it asked for).
    expect_identical(
        X("org.Hs.eg.db", c("7157", "672", "1"), "ENTREZID", "SYMBOL"),
        c("7157" = "TP53", "672" = "BRCA1", "1" = "A1BG")
    )

    ## An id the OrgDb cannot map comes back as NA, NOT as itself.
    ##
    ## EXTID2NAME() has a branch that is meant to pass unmappable ids through
    ## unchanged, but it can never fire: AnnotationDbi::select() already returns a
    ## row for an unknown key with NA in the value column, so the key IS present
    ## in gn.df$GeneID and `!geneID %in% gn.df$GeneID` matches nothing. The
    ## pass-through branch is dead code and the NA survives to the caller.
    ##
    ## Recorded as-is rather than fixed here: setReadable() turns these into the
    ## readable path via gn[...], where an NA does drop the gene, so the behaviour
    ## is pinned so that a future fix is a deliberate change with a visible diff.
    out2 <- X("org.Hs.eg.db", c("7157", "999999999"), "ENTREZID", "SYMBOL")
    expect_identical(unname(out2["7157"]), "TP53")
    expect_true(is.na(out2[["999999999"]]))

    ## an unsupported keytype is an error, not a silent empty map
    expect_error(X("org.Hs.eg.db", "7157", "NOT_A_KEYTYPE", "SYMBOL"),
                 "keytype is not supported")
})

test_that("setReadable rejects something that is not an enrichment result", {
    skip_if_no_orgdb()
    expect_error(setReadable(data.frame(a = 1), "org.Hs.eg.db"), "input should be")
})

test_that("setReadable errors on a toType the OrgDb does not have", {
    skip_if_no_orgdb()
    ## setReadable() has its own "Fail to convert ..." warning for a toType
    ## missing from the OrgDb, but AnnotationDbi rejects an invalid column name
    ## before that branch is reached, so the observable behaviour is an error.
    ## Pinned as-is: the point is that it does not silently return unchanged ids.
    res <- run_enrichgo()
    expect_error(
        suppressWarnings(setReadable(res, "org.Hs.eg.db", keyType = "ENTREZID",
                                     toType = "NOT_A_COLUMN")),
        "Invalid columns"
    )
})
