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
#' @param source One of `"tf"` (C3:TFT), `"mirna"` (C3:MIR),
#'   `"perturbagen"` (C2:CGP), `"cellstate"` (C8) or `"tissue"` (C4:CGN).
#' @param species Currently only human (`"Homo sapiens"`) is supported, because
#'   the files addressed here are the `*.Hs.*` release.
#' @param release MSigDB release string, used verbatim in the download URL and
#'   recorded in the `version` slot. Pinning it as a literal makes the
#'   collection contents reproducible; a release is never inferred from the
#'   current date.
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
mechanism_sets <- function(source = c("tf", "mirna", "perturbagen", "cellstate", "tissue"),
                           species = "Homo sapiens",
                           release = "2026.2.Hs",
                           cache_dir = getOption("clusterProfiler_cache_dir")) {
    source <- match.arg(source)

    if (!identical(species, "Homo sapiens")) {
        rlang::abort(
            c(
                paste0("Only `species = \"Homo sapiens\"` is supported, not \"", species, "\"."),
                i = "mechanism_sets() currently addresses the *.Hs.* MSigDB release only."
            )
        )
    }

    spec <- mechanism_set_files[[source]]
    file <- sub("<release>", release, spec$file, fixed = TRUE)
    cache_dir <- mechanism_cache_root(cache_dir)
    path <- file.path(cache_dir, file)

    accessed <- mechanism_fetch_gmt(path, release, cache_dir)
    parsed <- mechanism_parse_gmt(path)

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
