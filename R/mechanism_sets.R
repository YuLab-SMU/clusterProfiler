#' Upstream gene-set collections for mechanism analysis
#'
#' Download an MSigDB collection whose members are genes and whose set
#' identifiers name an *upstream participant* (a transcription factor, a
#' microRNA, a perturbagen, a cell state or a tissue) rather than a biological
#' process. Feeding such a collection to [enricher()] or [GSEA()] turns a gene
#' list into a rank of candidate mechanisms, which is the opposite direction
#' from a pathway analysis.
#'
#' The returned object is a [gson][gson::gson]-style `GSON`, accepted by
#' `enricher(..., gson = )` and `GSEA(..., gson = )`. It must go through the
#' `gson =` argument, not `TERM2GENE =`: only the `gson =` path carries
#' `species`, `keytype` and `gsname` into the result object, and
#' [evidence()] relies on those three slots.
#'
#' @param source One of `"tf"` (C3:TFT), `"tfunion"` (C3:TFT plus four public
#'   transcription-factor collections, see [mechanism_tfunion()]), `"mirna"` (C3:MIR),
#'   `"perturbagen"` (C2:CGP), `"cellstate"` (C8) or `"tissue"` (C4:CGN).
#' @param species Currently only human (`"Homo sapiens"`) is supported, because
#'   the files addressed here are the `*.Hs.*` release.
#' @param release MSigDB release string, used verbatim in the download URL and
#'   recorded in the `version` slot. Pinning it as a literal makes the
#'   collection contents reproducible; a release is never inferred from the
#'   current date.
#' @param subset Keep only the parts of the collection produced by a named
#'   evidence rule: `"gtrd"` for the ChIP-defined TF sets (whose identifiers are
#'   gene symbols), `"motif"` for the predicted ones, `"unnamed_motif"` for
#'   consensus sequences that name no factor. Measured on the airway list,
#'   `subset = "gtrd"` takes C3:TFT's share of family-or-consensus labels from
#'   40% to 1% -- but also its significant findings from 200 to 0, because the
#'   ChIP target sets overlap this list less than chance. Read before setting
#'   it.
#' @param cache_dir Directory holding downloaded GMT files. Defaults to
#'   `getOption("clusterProfiler_cache_dir")`, which [clusterProfiler-package]
#'   sets once per session.
#' @return A `GSON` object whose `gsid2gene` maps set identifiers to HGNC
#'   symbols, with `species`, `keytype = "SYMBOL"`, `gsname`, `version` and
#'   `accessed_date` filled in.
#' @seealso [evidence()] for the analysis that consumes these collections.
#'   `get_ppi_network()` for the network-side evidence source.
#' @examples
#' \dontrun{
#' tf <- mechanism_sets("tf")
#' x <- enricher(c("HIF1A", "PGK1", "CA9"), gson = tf)
#' }
#' @export
mechanism_sets <- function(source = c("tf", "tfunion", "mirna", "perturbagen", "cellstate", "tissue"),
                           species = "Homo sapiens",
                           release = "2026.2.Hs",
                           cache_dir = getOption("clusterProfiler_cache_dir"),
                           subset = NULL) {
    source <- match.arg(source)
    subset <- .mechanism_subset_source(source, subset)

    if (!identical(species, "Homo sapiens")) {
        rlang::abort(
            c(
                paste0("Only `species = \"Homo sapiens\"` is supported, not \"", species, "\"."),
                i = "mechanism_sets() currently addresses the *.Hs.* MSigDB release only."
            )
        )
    }

    if (identical(source, "tfunion")) {
        return(mechanism_tfunion(species = species, release = release,
                                 cache_dir = cache_dir))
    }
    spec <- mechanism_set_files[[source]]
    file <- sub("<release>", release, spec$file, fixed = TRUE)
    cache_dir <- mechanism_cache_root(cache_dir)
    path <- file.path(cache_dir, file)

    accessed <- mechanism_fetch_gmt(path, release, cache_dir)
    parsed <- mechanism_parse_gmt(path)

    ## A collection can hold two sources of very different quality under one
    ## name. C3:TFT is the case that matters: GTRD holds ChIP-defined targets
    ## for one gene each, while TFT_LEGACY holds motif matches, 73% of which
    ## name a family or a consensus sequence rather than a gene. Mixing them
    ## made 40% of the entities unqueryable.
    if (!is.null(subset)) {
        keep <- parsed$gsid2name$gsid[
            mechanism_rule_for(source, parsed$gsid2name$gsid) %in% subset]
        parsed$gsid2gene <- parsed$gsid2gene[parsed$gsid2gene$gsid %in% keep, ,
                                             drop = FALSE]
        parsed$gsid2name <- parsed$gsid2name[parsed$gsid2name$gsid %in% keep, ,
                                             drop = FALSE]
    }

    gson::gson(
        gsid2gene     = parsed$gsid2gene,
        gsid2name     = parsed$gsid2name,
        species       = species,
        gsname        = spec$gsname,
        version       = release,
        accessed_date = accessed,
        keytype       = "SYMBOL"
    )
}

## One row per collection: the file name pattern, and the gsname that
## enricher()/GSEA() will write into @ontology.
##
## The prefixes are literal because they cannot be derived from the collection
## name: c8 uses `c8.all`, c4:Cgn uses `c4.cgn`, and `c2.cgn` (the naming a
## reader would guess for a "C2 CGN" collection) does not exist -- MSigDB
## returns 404 for it. Verified against the 2026.2.Hs release.
## `complex` has no GMT file and is deliberately absent; it needs a hand-built
## GSON, so it is not offered here rather than silently missing.
mechanism_set_files <- list(
    tf = list(
        file   = "c3.tft.v<release>.symbols.gmt",
        gsname = "MSigDB C3:TFT"
    ),
    mirna = list(
        file   = "c3.mir.v<release>.symbols.gmt",
        gsname = "MSigDB C3:MIR"
    ),
    perturbagen = list(
        file   = "c2.cgp.v<release>.symbols.gmt",
        gsname = "MSigDB C2:CGP"
    ),
    cellstate = list(
        file   = "c8.all.v<release>.symbols.gmt",
        gsname = "MSigDB C8"
    ),
    tissue = list(
        file   = "c4.cgn.v<release>.symbols.gmt",
        gsname = "MSigDB C4:CGN"
    )
)

mechanism_gmt_url <- function(file, release) {
    paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",
           release, "/", file)
}

## Download `path` if it is absent or was fetched for a different release,
## and return the date it was first obtained.
##
## The release check matters: reusing a cached file after `release` changes
## would run the analysis against a collection the caller did not ask for,
## and the only symptom would be numbers that drift between sessions.
## The directory GMT files live in, resolved once.
##
## A test that wants to check against the real MSigDB download has to ask for
## the same place `mechanism_sets()` reads from; hard-coding a path worked only
## while the test was run from one particular working directory, which is not
## how R CMD check runs.
mechanism_cache_root <- function(cache_dir = NULL) {
    cache_dir %||% getOption("clusterProfiler_cache_dir") %||%
        yulab.utils::user_dir("clusterProfiler")
}

mechanism_fetch_gmt <- function(path, release, cache_dir) {
    meta <- paste0(path, ".meta")

    if (file.exists(path) && file.exists(meta)) {
        cached <- mechanism_read_meta(meta, "release")
        if (identical(cached, release)) {
            return(mechanism_read_meta(meta, "accessed_date"))
        }
        rlang::warn(paste0(
            "Cached ", basename(path), " was downloaded for release ",
            cached %||% "<unknown>", ", not ", release, "; downloading again."
        ))
    }

    url <- mechanism_gmt_url(basename(path), release)
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    status <- tryCatch(
        suppressWarnings(utils::download.file(url, path, quiet = TRUE)),
        error = function(e) conditionMessage(e)
    )
    if (!is.numeric(status) || !identical(as.integer(status), 0L) || !file.exists(path)) {
        rlang::abort(c(
            paste0("Could not download ", url),
            x = if (is.numeric(status)) paste0("download.file() returned ", status) else status,
            i = paste0("Destination: ", path)
        ))
    }
    ## A captive portal or proxy error page also arrives as a successful
    ## download, and would then parse into an empty collection.
    if (mechanism_is_html(path)) {
        unlink(path)
        rlang::abort(c(
            paste0("Downloaded an HTML page instead of a GMT file: ", url),
            i = "A proxy or login page is intercepting the request."
        ))
    }

    accessed <- as.character(Sys.Date())
    writeLines(c(paste0("release=", release), paste0("accessed_date=", accessed)), meta)
    accessed
}

mechanism_is_html <- function(path) {
    bytes <- readBin(path, "raw", n = 512L)
    bytes <- bytes[!bytes %in% charToRaw(" \t\r\n")]
    length(bytes) > 0 && bytes[1] == charToRaw("<")
}

mechanism_read_meta <- function(meta, key) {
    lines <- readLines(meta, warn = FALSE)
    hit <- lines[startsWith(lines, paste0(key, "="))]
    if (!length(hit)) return(NULL)
    sub(paste0("^", key, "="), "", hit[1])
}

## Parse a GMT file into the two tables gson::gson() takes.
##
## GMT layout: set name <TAB> source URL <TAB> member 1 <TAB> member 2 ...
## Line endings are stripped explicitly because MSigDB files have been
## distributed with CRLF.
mechanism_parse_gmt <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- sub("\r$", "", lines)
    lines <- lines[nzchar(lines)]
    parts <- strsplit(lines, "\t", fixed = TRUE)

    gsid2gene <- do.call(rbind, lapply(parts, function(x) {
        data.frame(gsid = x[1], gene = x[-c(1, 2)], stringsAsFactors = FALSE)
    }))
    gsid2gene <- unique(gsid2gene[!is.na(gsid2gene$gene), , drop = FALSE])

    if (nrow(gsid2gene) == 0) {
        rlang::abort(paste0("No gene sets could be parsed from ", path))
    }

    gsid2name <- unique(data.frame(
        gsid = gsid2gene$gsid,
        name = gsid2gene$gsid,
        stringsAsFactors = FALSE
    ))

    list(gsid2gene = gsid2gene, gsid2name = gsid2name)
}

## Restrict a collection to the parts worth using.
##
## `subset` is expressed in the vocabulary of the xref rules rather than in
## set-name patterns, so a caller can say "only what was measured" and not have
## to know that GTRD happens to be named `*_TARGET_GENES`.
.mechanism_subset_source <- function(source, subset) {
    if (is.null(subset)) return(NULL)
    allowed <- c("gtrd", "motif", "unnamed_motif", "predicted_target", "signature",
                 "recorded_perturbation")
    bad <- setdiff(subset, allowed)
    if (length(bad)) {
        rlang::abort(c(
            paste0("Unknown `subset`: ", paste(bad, collapse = ", "), "."),
            i = paste0("Known rules: ", paste(allowed, collapse = ", "), ".")
        ))
    }
    subset
}

#' The transcription-factor union collection
#'
#' @description
#' Builds one collection from C3:TFT plus four public transcription-factor
#' libraries, so the search for an upstream regulator is not bounded by the
#' ~592 factors C3:TFT happens to name.
#'
#' @details
#' **Why this exists.** C3:TFT has no target set at all for `EPAS1`, `NFKB1`,
#' `TP53`, `CTCF`, `ESR1`, `NOTCH1`, `JUN` or `FOS`. On a hypoxia dataset whose
#' published answer is HIF-2 (`EPAS1`), no amount of grading can find a
#' regulator the collection does not name -- the ceiling was the data, not the
#' method. These four libraries together add 361 factors C3:TFT lacks.
#'
#' **Why the sets stay separate.** Merging each factor's sets into one gene list
#' would be wrong, and measurably so: CTCF's 108 ENCODE sets union to 21295
#' genes, 11% of the background, which no enrichment test can interpret. Each
#' set keeps its identity and is tagged with its library prefix
#' (`ENCODE:`, `ChEA:`, `TFP:`, `TRRUST:`, `C3TFT:`); folding duplicates is
#' [consensus()]'s job, and it has the evidence levels to do that with.
#'
#' **What it costs.** 3258 sets against C3:TFT's 1114, so the multiple-testing
#' correction is spread over three times as many hypotheses. Measured on the
#' hypoxia contrast, the true signal survived that comfortably -- `HIF1A` moved
#' from padj 1.2e-4 (C3:TFT) to 9.3e-28 -- while three random gene sets of the
#' same size returned **zero** significant sets, so the extra hypotheses did not
#' buy extra false positives.
#'
#' @inheritParams mechanism_sets
#' @return A `GSON`, in the same shape as [mechanism_sets()] returns.
#' @seealso [mechanism_sets()]
#' @examples
#' \dontrun{
#' ev <- evidence(genes, upstream = "tfunion")
#' }
#' @export
mechanism_tfunion <- function(species = "Homo sapiens",
                              release = "2026.2.Hs",
                              cache_dir = getOption("clusterProfiler_cache_dir")) {
    if (!identical(species, "Homo sapiens")) {
        rlang::abort("The TF union collection is human-only: the libraries are named for Homo sapiens.")
    }
    cache_dir <- mechanism_cache_root(cache_dir)

    base <- mechanism_sets("tf", species = species, release = release,
                           cache_dir = cache_dir)
    parts <- list(data.frame(gsid = paste0("C3TFT:", base@gsid2gene$gsid),
                             gene = base@gsid2gene$gene,
                             stringsAsFactors = FALSE))
    gsid2name <- data.frame(gsid = paste0("C3TFT:", base@gsid2name$gsid),
                            name = base@gsid2name$name,
                            stringsAsFactors = FALSE)

    for (tag in names(mechanism_tfunion_parts)) {
        spec <- mechanism_tfunion_parts[[tag]]
        df <- tryCatch(.mechanism_enrichr_gmt(spec$library, cache_dir, tag),
                       error = function(e) {
                           rlang::warn(paste0("TF union: could not fetch ",
                                              spec$library, " (",
                                              conditionMessage(e), "). The union ",
                                              "proceeds without it."))
                           NULL
                       })
        if (is.null(df)) next
        parts[[length(parts) + 1L]] <- df
        gsid2name <- rbind(gsid2name, data.frame(
            gsid = unique(df$gsid),
            name = sub("^[^:]+:", "", unique(df$gsid)),
            stringsAsFactors = FALSE))
    }

    gsid2gene <- unique(do.call(rbind, parts))
    ## Ensembl ids and retired symbols are in some of these libraries; a gene
    ## that cannot be matched to a symbol cannot overlap a caller's list, and
    ## leaving it in inflates setSize and dilutes every ratio.
    known <- mechanism_known_symbols(species)
    gsid2gene <- gsid2gene[gsid2gene$gene %in% known, , drop = FALSE]

    gson::gson(
        gsid2gene     = gsid2gene,
        gsid2name     = gsid2name[!duplicated(gsid2name$gsid), , drop = FALSE],
        species       = species,
        gsname        = "TF union (C3:TFT + ENCODE + ChEA + TF_Perturbations + TRRUST)",
        version       = release,
        accessed_date = base@accessed_date,
        keytype       = "SYMBOL"
    )
}

## Fetch one Enrichr library as a GSON-shaped data frame.
##
## Enrichr serves plain GMT with a blank second column, which is why the gene
## list drops two fields rather than one.
.mechanism_enrichr_gmt <- function(library, cache_dir, tag,
                                   timeout = getOption("clusterProfiler_timeout", 120)) {
    file <- paste0("enrichr.", gsub("[^A-Za-z0-9_]", "_", library), ".gmt")
    path <- file.path(cache_dir, file)
    if (!file.exists(path)) {
        url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
                      utils::URLencode(library, reserved = TRUE))
        dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
        status <- tryCatch(
            suppressWarnings(utils::download.file(url, path, quiet = TRUE)),
            error = function(e) conditionMessage(e))
        if (!is.numeric(status) || !identical(as.integer(status), 0L) || !file.exists(path)) {
            rlang::abort(paste0("Could not download ", library))
        }
        if (mechanism_is_html(path)) {
            unlink(path)
            rlang::abort(paste0("Got an HTML page instead of a GMT for ", library))
        }
    }
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- lines[nzchar(lines)]
    if (!length(lines)) {
        rlang::abort(paste0("The library file for ", library, " is empty."))
    }
    sp <- strsplit(lines, "\t", fixed = TRUE)
    gsid <- vapply(sp, function(x) x[1], "")
    genes <- lapply(sp, function(x) unique(x[-(1:2)]))
    ## Every one of these libraries mixes species inside the set identifier,
    ## and the mix is not small: 49% of ChEA's sets are Mouse or Rat and 17% of
    ## ENCODE's are mm9. A mouse target list matched against a human gene set
    ## produces a confident wrong answer, so the species is decided here rather
    ## than left to the symbol overlap to filter out.
    keep <- switch(tag,
        TRRUST = grepl(" human$", gsid),
        TFP    = grepl("HUMAN", gsid) & grepl(" (KD|KO) ", gsid),
        ChEA   = grepl(" Human", gsid) & !grepl(" Mouse| Rat", gsid),
        ENCODE = grepl("hg19|hg38", gsid) & !grepl("mm9|mm10|mm39|rat", gsid),
        rep(TRUE, length(gsid)))
    gsid <- gsid[keep]; genes <- genes[keep]
    ids <- paste0(tag, ":", gsid)
    data.frame(gsid = rep(ids, lengths(genes)),
               gene = unlist(genes, use.names = FALSE),
               stringsAsFactors = FALSE)
}

## The transcription-factor union: which library each prefix is, and what its
## evidence is worth.
##
## The kinds differ deliberately. ENCODE and ChEA assign genes to ChIP peaks,
## which measures binding. TFP records which genes moved when a factor was
## knocked down, which measures an effect but not a direct target. TRRUST is a
## literature curation with no assay visible to us, so it is reported as an
## observation and cannot corroborate anything on its own.
mechanism_tfunion_parts <- list(
    ENCODE = list(library = "ENCODE_TF_ChIP-seq_2015",
                  rule = "chip_union", level = "B"),
    ChEA   = list(library = "ChEA_2022", rule = "chip_union", level = "B"),
    TFP    = list(library = "TF_Perturbations_Followed_by_Expression",
                  rule = "perturb_signature", level = "B"),
    TRRUST = list(library = "TRRUST_Transcription_Factors_2019",
                  rule = "curated_target", level = "C")
)

## Which library prefix a union set identifier carries, or NA for a plain id.
mechanism_union_prefix <- function(gsid) {
    ## `[A-Za-z0-9]`, not `[A-Z0-9]`: the ChEA prefix carries a lowercase
    ## letter. An uppercase-only pattern does not error on it -- it routes those
    ## 653 sets down the C3:TFT naming path and silently mangles their names.
    hit <- grepl("^[A-Za-z0-9]+:", gsid)
    out <- rep(NA_character_, length(gsid))
    ## strsplit rather than a regex back-reference: a capture written as a
    ## backslash-one in R source needs two backslashes, and one short yields the
    ## control character 0x01 instead of the capture. That is not an error -- it
    ## is a column of unprintable strings that then fails every symbol match.
    out[hit] <- vapply(strsplit(gsid[hit], ":", fixed = TRUE),
                       function(x) x[1], "", USE.NAMES = FALSE)
    out
}

## The regulator a union set identifier names.
##
## The four added libraries all put the factor first and follow it with cell
## line, assay or accession, so the first whitespace-delimited token is the
## regulator. C3:TFT sets keep their own, harder naming (motif prefixes and
## quality tiers), so they are routed to mechanism_regulator() instead.
mechanism_union_entity <- function(gsid, species = "Homo sapiens") {
    gsid <- as.character(gsid)
    prefix <- mechanism_union_prefix(gsid)
    bare <- sub("^[^:]+:", "", gsid)
    out <- character(length(gsid))
    is_c3 <- !is.na(prefix) & prefix == "C3TFT"
    is_other <- !is.na(prefix) & !is_c3
    plain <- is.na(prefix)
    if (any(is_c3)) out[is_c3] <- mechanism_regulator(bare[is_c3],
                                                      symbols = mechanism_known_symbols(species))
    if (any(plain)) out[plain] <- mechanism_regulator(bare[plain],
                                                      symbols = mechanism_known_symbols(species))
    if (any(is_other)) {
        tok <- vapply(strsplit(bare[is_other], " ", fixed = TRUE),
                         function(x) x[1], "", USE.NAMES = FALSE)
        ## ENCODE writes the phospho form of POLR2A as one run-on word.
        out[is_other] <- sub("phosphoS2$", "", tok)
    }
    out
}
