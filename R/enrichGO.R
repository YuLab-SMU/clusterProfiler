#' GO Enrichment Analysis of a gene set.
#' Given a vector of genes, this function will return the enrichment GO
#' categories after FDR control.
#'
#'
#' @param gene a vector of entrez gene id.
#' @param OrgDb OrgDb
#' @param keyType keytype of input gene
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @inheritParams enricher
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param readable whether mapping gene ID to gene Name
#' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
#' @param evidence optional character vector of GO evidence codes to keep (for
#'   example `c("IDA", "IPI", "IMP")`, or `"IEA"` for electronic annotations
#'   only). When supplied, annotations carrying any other evidence code are
#'   dropped before the analysis — useful to avoid circular reasoning (dropping
#'   `IEP` when the gene sets come from co-expression clustering) or to keep only
#'   manually reviewed annotations (dropping `IEA`). Because the annotation
#'   changes, so does the background, so the resulting terms are not simply a
#'   subset of the unfiltered run. Filtered runs are not cached, so they rebuild
#'   the annotation each time. `NULL` (default) uses every annotation.
#' @return An \code{enrichResult} instance.
#' @importClassesFrom enrichit enrichResult
#' @importFrom enrichit setReadable
#' @importFrom utils stack
#' @seealso [enrichResult-class], [compareCluster]
#' @keywords manip
#' @export
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @examples
#' \dontrun{
#'   data(geneList, package = "DOSE")
#' 	de <- names(geneList)[1:100]
#' 	yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
#' 	head(yy)
#' }
enrichGO <- function(
    gene,
    OrgDb,
    keyType = "ENTREZID",
    ont = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE,
    evidence = NULL
) {
    has_universe <- !missing(universe)
    universe_value <- if (has_universe) universe else NULL

    # This is to avoid Memory Boom for non-ENTREZID keyType
    # see https://github.com/YuLab-SMU/clusterProfiler/issues/805
    #
    # If we use `get_GO_data` to build the GSON object for non-ENTREZID keyType,
    # it will trigger `AnnotationDbi::mapIds` to map all GO terms to the specific keyType.
    # This will be extremely slow and memory consuming if the keyType is ACCNUM,
    # because one gene can have multiple ACCNUMs.
    #
    # So if the OrgDb supports ENTREZID, we first map the input gene to ENTREZID,
    # and then use ENTREZID to do the enrichment analysis.
    # After that, we map the result back to the original keyType.
    #
    # If the OrgDb does NOT have ENTREZID (e.g., some custom OrgDb for non-model
    # organisms), we skip this optimization and use the original keyType directly.
    # see https://github.com/YuLab-SMU/clusterProfiler/issues/823

    ont <- toupper(ont)
    ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

    if (keyType != "ENTREZID" &&
        "ENTREZID" %in% keytypes(yulab.utils::load_OrgDb(OrgDb))) {
        OrgDb <- yulab.utils::load_OrgDb(OrgDb)
        gene_ids <- unique(as.character(gene))
        gene_map <- map_keytype_to_entrez(gene_ids, keyType, OrgDb)
        gene_complete <- all(gene_ids %in% gene_map$source)

        universe_map <- NULL
        universe_complete <- TRUE
        if (has_universe && !is.null(universe)) {
            universe_ids <- unique(as.character(universe))
            universe_map <- map_keytype_to_entrez(universe_ids, keyType, OrgDb)
            universe_complete <- all(universe_ids %in% universe_map$source)
        }

        # Keep the #805 fast path whenever every requested ID maps.  The
        # hybrid path below is only needed when a source ID has no ENTREZID.
        if (gene_complete && universe_complete) {
            gene_entrez <- unique(gene_map$entrez[gene_map$source %in% gene_ids])
            if (!length(gene_entrez)) {
                message("--> No gene can be mapped....")
                return(NULL)
            }

            universe_entrez <- NULL
            if (has_universe && !is.null(universe)) {
                universe_entrez <- unique(
                    universe_map$entrez[universe_map$source %in% universe_ids]
                )
                if (!length(universe_entrez)) {
                    message("--> No universe can be mapped....")
                    return(NULL)
                }
            }

            res <- enrichGO(
                gene = gene_entrez,
                OrgDb = OrgDb,
                keyType = "ENTREZID",
                ont = ont,
                pvalueCutoff = pvalueCutoff,
                pAdjustMethod = pAdjustMethod,
                universe = universe_entrez,
                qvalueCutoff = qvalueCutoff,
                minGSSize = minGSSize,
                maxGSSize = maxGSSize,
                readable = readable,
                pool = pool,
                evidence = evidence
            )

            if (is.null(res)) {
                return(NULL)
            }

            # setReadable() writes @keytype itself.  Restore the user's
            # original key type after the conversion (#829).
            if (readable) {
                res <- setReadable(res, OrgDb, keyType = "ENTREZID")
            } else {
                res <- setReadable(
                    res,
                    OrgDb,
                    keyType = "ENTREZID",
                    toType = keyType
                )
                res@readable <- FALSE
            }
            res@keytype <- keyType
            res <- restore_entrez_result_ids(
                res,
                OrgDb = OrgDb,
                keyType = keyType,
                original_gene = gene_ids,
                original_universe = if (has_universe) universe_ids else NULL
            )
            return(res)
        }

        # At least one source ID cannot be represented by ENTREZID.  Use one
        # unified GSON containing the Entrez annotation plus direct native
        # annotations for the source-only IDs, rather than rebuilding the full
        # native GO -> keyType table (which is the memory-heavy #805 path).
        if (ont == "ALL" && !pool) {
            lres <- lapply(c("BP", "CC", "MF"), function(one_ont) {
                suppressMessages(enrichGO(
                    gene = gene,
                    OrgDb = OrgDb,
                    keyType = keyType,
                    ont = one_ont,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod = pAdjustMethod,
                    universe = universe_value,
                    qvalueCutoff = qvalueCutoff,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    readable = readable,
                    pool = TRUE,
                    evidence = evidence
                ))
            })
            lres <- lres[!vapply(lres, is.null, logical(1))]
            if (!length(lres)) {
                return(NULL)
            }
            df <- do.call("rbind", lapply(lres, as.data.frame))
            geneSets <- lres[[1]]@geneSets
            if (length(lres) > 1) {
                for (i in 2:length(lres)) {
                    geneSets <- append(geneSets, lres[[i]]@geneSets)
                }
            }
            res <- lres[[1]]
            res@result <- df
            res@geneSets <- geneSets
            return(res)
        }

        hybrid <- build_hybrid_go_data(
            OrgDb = OrgDb,
            ont = ont,
            keyType = keyType,
            gene = gene_ids,
            universe = if (has_universe) universe else NULL,
            evidence = evidence
        )

        gene_internal <- hybrid_ids_for_input(gene_ids, hybrid$mapping, keyType)
        universe_internal <- NULL
        if (has_universe && !is.null(universe)) {
            universe_internal <- hybrid_ids_for_input(
                unique(as.character(universe)), hybrid$mapping, keyType
            )
        } else {
            universe_internal <- unique(hybrid$gson@gsid2gene$gene)
            message(
                "--> No universe supplied; using ",
                length(universe_internal),
                " genes with ", ont,
                " annotations as the background (hybrid Entrez/native annotation)."
            )
        }

        res <- enrichit::ora_gson(
            gene_internal,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod,
            universe = universe_internal,
            qvalueCutoff = qvalueCutoff,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            gson = hybrid$gson
        )
        if (is.null(res)) {
            return(NULL)
        }

        res@organism <- get_organism(OrgDb)
        res@ontology <- ont
        res@keytype <- keyType
        res <- restore_hybrid_result(
            res,
            OrgDb = OrgDb,
            keyType = keyType,
            mapping = hybrid$mapping,
            readable = readable,
            original_gene = gene_ids,
            original_universe = if (has_universe) universe else NULL
        )
        if (ont == "ALL") {
            res <- add_GO_Ontology(res, hybrid$gson)
        }
        return(res)
    }

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))
    GO_DATA <- get_GO_data(OrgDb, ont, keyType, evidence = evidence)

    if (missing(universe)) {
        universe <- NULL
    }

    if (ont == "ALL" && !pool) {
        lres <- lapply(c("BP", "CC", "MF"), function(ont) {
            suppressMessages(enrichGO(
                gene,
                OrgDb,
                keyType,
                ont,
                pvalueCutoff,
                pAdjustMethod,
                universe,
                qvalueCutoff,
                minGSSize,
                maxGSSize,
                evidence = evidence
            ))
        })

        lres <- lres[!vapply(lres, is.null, logical(1))]
        if (length(lres) == 0) {
            return(NULL)
        }

        df <- do.call('rbind', lapply(lres, as.data.frame))
        geneSets <- lres[[1]]@geneSets
        if (length(lres) > 1) {
            for (i in 2:length(lres)) {
                geneSets <- append(geneSets, lres[[i]]@geneSets)
            }
        }
        res <- lres[[1]]
        res@result <- df
        res@geneSets <- geneSets
    } else {
        res <- enrichit::ora_gson(
            gene,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod,
            universe = universe,
            qvalueCutoff = qvalueCutoff,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            gson = GO_DATA
        )

        if (is.null(res)) {
            return(res)
        }
    }
    if (keyType == 'SYMBOL') {
        res@readable <- TRUE
    }
    res@keytype <- keyType
    res@organism <- get_organism(OrgDb)
    if (readable) {
        res <- setReadable(res, OrgDb)
    }
    res@ontology <- ont

    if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA)
    }
    return(res)
}

map_keytype_to_entrez <- function(ids, fromType, OrgDb, chunk_size = 50000L) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids)]
    if (!length(ids)) {
        return(data.frame(source = character(), entrez = character(),
                          stringsAsFactors = FALSE))
    }

    chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))
    mapped <- lapply(chunks, function(chunk) {
        x <- suppressWarnings(suppressMessages(
            AnnotationDbi::select(
                OrgDb,
                keys = chunk,
                keytype = fromType,
                columns = "ENTREZID"
            )
        ))
        if (!nrow(x) || !fromType %in% names(x) || !"ENTREZID" %in% names(x)) {
            return(NULL)
        }
        x <- x[, c(fromType, "ENTREZID"), drop = FALSE]
        colnames(x) <- c("source", "entrez")
        x
    })
    mapped <- mapped[!vapply(mapped, is.null, logical(1))]
    if (!length(mapped)) {
        return(data.frame(source = character(), entrez = character(),
                          stringsAsFactors = FALSE))
    }

    mapped <- unique(do.call(rbind, mapped))
    mapped <- mapped[!is.na(mapped$source) & !is.na(mapped$entrez), , drop = FALSE]
    mapped$source <- as.character(mapped$source)
    mapped$entrez <- as.character(mapped$entrez)
    mapped
}

hybrid_ids_for_input <- function(ids, mapping, keyType) {
    ids <- unique(as.character(ids))
    if (!length(ids)) {
        return(character())
    }

    ## One split, then a named lookup. Comparing `mapping$source == id` per id
    ## is linear in the whole mapping for every input id; on a 20k-row mapping
    ## and ~600 ids per cluster that comparison alone was 7.8 s across 13
    ## clusters (2026-10-18). A source id may map to several Entrez ids, so the
    ## lookup keeps every match, in mapping order, as before.
    by_source <- split(mapping$entrez, mapping$source)
    out <- lapply(ids, function(id) {
        entrez <- unique(by_source[[id]])
        if (length(entrez)) {
            paste0("E:", entrez)
        } else {
            paste0("N:", keyType, ":", id)
        }
    })
    unique(unlist(out, use.names = FALSE))
}

build_hybrid_go_data <- function(OrgDb, ont, keyType, gene, universe = NULL,
                                 evidence = NULL) {
    OrgDb <- yulab.utils::load_OrgDb(OrgDb)
    org <- get_organism(OrgDb)
    ## Cache only the universe-free, unfiltered build: with a universe the
    ## source ids include the gene list and change per call, and an
    ## evidence-filtered annotation must never share a cache slot with the
    ## unfiltered one (same rule as get_GO_data()).
    cacheable <- is.null(universe) && is.null(evidence)
    GO_Env <- get_GO_Env()
    if (cacheable && hybrid_cache_usable(GO_Env, org, ont, keyType)) {
        return(get("hybrid_go", envir = GO_Env))
    }
    has_universe <- !is.null(universe)
    source_ids <- if (has_universe) {
        unique(c(as.character(gene), as.character(universe)))
    } else {
        AnnotationDbi::keys(OrgDb, keytype = keyType)
    }

    mapping <- map_keytype_to_entrez(source_ids, keyType, OrgDb)
    native_ids <- setdiff(source_ids, unique(mapping$source))

    entrez_gson <- get_GO_data(
        OrgDb,
        ont = ont,
        keytype = "ENTREZID",
        evidence = evidence
    )
    entrez_anno <- unique(entrez_gson@gsid2gene[, c("gsid", "gene"), drop = FALSE])
    entrez_anno$gene <- paste0("E:", as.character(entrez_anno$gene))

    native_anno <- data.frame(
        gsid = character(),
        gene = character(),
        stringsAsFactors = FALSE
    )
    if (length(native_ids)) {
        columns <- c("GOALL", "ONTOLOGYALL")
        if (!is.null(evidence)) {
            columns <- c(columns, "EVIDENCEALL")
        }
        native <- tryCatch(
            suppressWarnings(suppressMessages(
                AnnotationDbi::select(
                    OrgDb,
                    keys = native_ids,
                    keytype = keyType,
                    columns = columns
                )
            )),
            error = function(e) data.frame()
        )
        if (nrow(native)) {
            native <- native[!is.na(native$GOALL) &
                             !is.na(native$ONTOLOGYALL), , drop = FALSE]
            if (ont != "ALL") {
                native <- native[native$ONTOLOGYALL == ont, , drop = FALSE]
            }
            if (!is.null(evidence) && "EVIDENCEALL" %in% names(native)) {
                native <- native[native$EVIDENCEALL %in% evidence, , drop = FALSE]
            }
            if (nrow(native)) {
                native_anno <- unique(data.frame(
                    gsid = as.character(native$GOALL),
                    gene = paste0("N:", keyType, ":",
                                  as.character(native[[keyType]])),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    gsid2gene <- unique(rbind(entrez_anno, native_anno))
    if (!nrow(gsid2gene)) {
        stop("No GO annotation is available for the requested hybrid key space.")
    }

    gsid2name <- entrez_gson@gsid2name
    if (!is.null(gsid2name)) {
        gsid2name <- gsid2name[
            gsid2name$gsid %in% unique(gsid2gene$gsid),
            , drop = FALSE
        ]
    }

    hybrid_gson <- gson::gson(
        gsid2gene = gsid2gene,
        gsid2name = gsid2name,
        species = entrez_gson@species,
        gsname = entrez_gson@gsname,
        version = entrez_gson@version,
        accessed_date = entrez_gson@accessed_date,
        keytype = keyType
    )

    result <- list(gson = hybrid_gson, mapping = mapping)
    if (cacheable) {
        assign("hybrid_go", result, envir = GO_Env)
        assign("hybrid_organism", org, envir = GO_Env)
        assign("hybrid_keytype", keyType, envir = GO_Env)
        assign("hybrid_ont", ont, envir = GO_Env)
    }
    result
}

hybrid_source_labels <- function(ids, mapping, keyType) {
    ids <- as.character(ids)
    labels <- ids
    native <- startsWith(ids, paste0("N:", keyType, ":"))
    labels[native] <- sub(paste0("^N:", keyType, ":"), "", ids[native])

    entrez <- startsWith(ids, "E:")
    if (any(entrez)) {
        e <- sub("^E:", "", ids[entrez])
        reverse <- mapping$source[match(e, mapping$entrez)]
        reverse[is.na(reverse)] <- e[is.na(reverse)]
        labels[entrez] <- reverse
    }
    names(labels) <- ids
    labels
}

hybrid_symbol_labels <- function(ids, source_labels, mapping, OrgDb, keyType) {
    symbols <- source_labels[ids]
    entrez <- startsWith(ids, "E:")
    native <- startsWith(ids, paste0("N:", keyType, ":"))

    if (any(native)) {
        native_ids <- source_labels[ids[native]]
        native_names <- tryCatch(
            enrichit::EXTID2NAME(OrgDb, native_ids, keyType, "SYMBOL"),
            error = function(e) setNames(native_ids, native_ids)
        )
        symbols[native] <- unname(native_names[native_ids])
    }

    if (any(entrez)) {
        e <- sub("^E:", "", ids[entrez])
        mapped_source <- mapping$source[match(e, mapping$entrez)]
        have_source <- !is.na(mapped_source)
        if (any(have_source)) {
            source_names <- tryCatch(
                enrichit::EXTID2NAME(
                    OrgDb,
                    mapped_source[have_source],
                    keyType,
                    "SYMBOL"
                ),
                error = function(e) setNames(
                    mapped_source[have_source],
                    mapped_source[have_source]
                )
            )
            symbols[which(entrez)[have_source]] <-
                unname(source_names[mapped_source[have_source]])
        }
        if (any(!have_source)) {
            entrez_names <- tryCatch(
                enrichit::EXTID2NAME(
                    OrgDb,
                    e[!have_source],
                    "ENTREZID",
                    "SYMBOL"
                ),
                error = function(err) setNames(e[!have_source], e[!have_source])
            )
            symbols[which(entrez)[!have_source]] <-
                unname(entrez_names[e[!have_source]])
        }
    }
    names(symbols) <- ids
    symbols
}

restore_entrez_result_ids <- function(x, OrgDb, keyType, original_gene,
                                      original_universe = NULL) {
    internal <- unique(unlist(x@geneSets, use.names = FALSE))
    internal <- as.character(internal[!is.na(internal)])
    source_labels <- if (length(internal)) {
        names_map <- tryCatch(
            enrichit::EXTID2NAME(OrgDb, internal, "ENTREZID", keyType),
            error = function(e) setNames(internal, internal)
        )
        setNames(unname(names_map), internal)
    } else {
        character()
    }

    x@geneSets <- map_gene_set_labels(x@geneSets, source_labels)$sets
    result_gene_ids <- unique(unlist(strsplit(
        as.character(x@result$geneID), "/", fixed = TRUE
    )))
    x@gene <- unique(c(
        as.character(original_gene),
        unlist(x@geneSets, use.names = FALSE),
        result_gene_ids
    ))

    if (!is.null(original_universe)) {
        x@universe <- unique(as.character(original_universe))
    } else if (length(x@universe)) {
        universe_ids <- as.character(x@universe)
        universe_map <- tryCatch(
            enrichit::EXTID2NAME(
                OrgDb,
                universe_ids,
                "ENTREZID",
                keyType
            ),
            error = function(e) setNames(universe_ids, universe_ids)
        )
        x@universe <- unique(unname(universe_map))
    }
    x@keytype <- keyType
    x
}
restore_hybrid_result <- function(x, OrgDb, keyType, mapping, readable,
                                  original_gene, original_universe = NULL) {
    internal <- unique(c(
        as.character(x@gene),
        as.character(x@universe),
        unlist(x@geneSets, use.names = FALSE)
    ))
    internal <- internal[!is.na(internal)]
    source_labels <- hybrid_source_labels(internal, mapping, keyType)
    symbol_labels <- hybrid_symbol_labels(
        internal, source_labels, mapping, OrgDb, keyType
    )

    gc <- geneInCategory(x)
    if (length(gc)) {
        values <- if (readable || keyType == "SYMBOL") {
            symbol_labels
        } else {
            source_labels
        }
        x@result$geneID <- map_gene_set_labels(gc, values)$gene_id
    }

    x@geneSets <- map_gene_set_labels(x@geneSets, source_labels)$sets
    result_gene_ids <- unique(unlist(strsplit(
        as.character(x@result$geneID), "/", fixed = TRUE
    )))
    x@gene <- unique(c(
        unname(source_labels[as.character(x@gene)]),
        unlist(x@geneSets, use.names = FALSE),
        result_gene_ids
    ))
    x@universe <- unique(unname(source_labels[as.character(x@universe)]))
    x@gene2Symbol <- symbol_labels
    x@keytype <- keyType
    x@readable <- isTRUE(readable) || keyType == "SYMBOL"
    x
}

#' Map the gene IDs of many gene sets in one pass
#'
#' `restore_hybrid_result()` and `restore_entrez_result_ids()` both translate
#' the IDs inside every gene set back to the source keyType. Doing that with
#' one `labels[ids]` + `unname()` per set is cheap per set and expensive in
#' aggregate: on a BP result with ~5,900 tested sets (390k ID slots) it cost
#' 7.5 s per call, measured 2026-10-18, and `compareCluster()` calls it once per
#' cluster. This maps the concatenated IDs once, splits the result back by the
#' original set lengths, and preserves both the list names (GO ids) and the
#' empty sets. Output is byte-identical to the per-set loop (pinned by a test).
#'
#' @param ids_list named list of character ID vectors, one per gene set
#' @param labels named character vector mapping internal ID -> label
#' @return list(sets = named list of unique labels, gene_id = named vector of
#'   "/"-joined labels)
#' @noRd
map_gene_set_labels <- function(ids_list, labels) {
    if (!length(ids_list)) {
        return(list(sets = list(), gene_id = character()))
    }
    all_ids <- unlist(ids_list, use.names = FALSE)
    mapped <- unname(labels[as.character(all_ids)])
    group <- rep(seq_along(ids_list), lengths(ids_list))
    parts <- split(mapped, factor(group, levels = seq_along(ids_list)), drop = FALSE)
    sets <- lapply(parts, unique)
    names(sets) <- names(ids_list)
    gene_id <- vapply(parts, function(ids) {
        paste0(unique(ids), collapse = "/")
    }, character(1))
    names(gene_id) <- names(ids_list)
    list(sets = sets, gene_id = gene_id)
}

#' Can the cached GO annotation serve this request?
#'
#' Compares a requested `(OrgDb organism, keytype)` against what is already
#' in the package's GO data environment.
#'
#' The cache always holds the **all-ontology** annotation, and `ont` is
#' deliberately not part of the key: a stored all-ontology annotation is a
#' superset of every single-ontology view, which is filtered per call by
#' `get_GO_data()`. Keying the cache on `ont` evicted the annotation whenever a
#' session switched ontology — BP then MF then BP paid the full
#' `AnnotationDbi::mapIds(GOALL → keytype)` rebuild twice, 22.3 s per miss on
#' human SYMBOL (2026-10-27); five alternating BP/MF `enrichGO()` calls took
#' 85.5 s and 11.9 s after the key change.
#'
#' Every comparison goes through `isTRUE()`: `get_organism()` returns `NA` for an
#' OrgDb it cannot identify — an AnnotationHub OrgDb, for instance — and a bare
#' `if (NA && ...)` aborts the whole call with the unhelpful
#' "missing value where TRUE/FALSE needed" (#232). Treating an unidentifiable
#' organism as "not cached" is the correct fallback: the annotation is rebuilt.
#'
#' @param GO_Env the package GO data environment
#' @param org organism of the requested OrgDb, possibly NA
#' @param ont requested ontology; accepted for call compatibility, not compared
#' @param keytype requested keytype
#' @return TRUE when the cache can be reused
#' @noRd
go_cache_usable <- function(GO_Env, org, ont, keytype) {
    for (what in c("organism", "keytype", "goAnno")) {
        if (!exists(what, envir = GO_Env, inherits = FALSE)) {
            return(FALSE)
        }
    }

    isTRUE(org == get("organism", envir = GO_Env)) &&
        isTRUE(keytype == get("keytype", envir = GO_Env))
}

#' Can a cached hybrid GO data object be reused?
#'
#' `build_hybrid_go_data()` maps every OrgDb key to Entrez and merges the
#' Entrez and native (unmapped) annotations. When no `universe` is supplied the
#' input is `AnnotationDbi::keys(OrgDb, keytype)` -- independent of the gene
#' list -- so the result depends only on (organism, ont, keytype). Without this
#' cache the whole build is repeated once per cluster: 7.7 s per call on human
#' GO BP (PBMC 10k, 2026-10-18), with `restore_hybrid_result()` adding another
#' 7.2 s.
#'
#' Deliberately stricter than `go_cache_usable()`: the ontology must match
#' exactly. The Entrez annotation cache may hold `ont = "ALL"` and be filtered
#' afterwards, but the hybrid object is already built and has no later filter.
#'
#' @param GO_Env the package GO data environment
#' @param org organism of the requested OrgDb, possibly NA
#' @param ont requested ontology
#' @param keytype requested keytype
#' @return TRUE when the cache can be reused
#' @noRd
hybrid_cache_usable <- function(GO_Env, org, ont, keytype) {
    for (what in c("hybrid_go", "hybrid_organism", "hybrid_keytype", "hybrid_ont")) {
        if (!exists(what, envir = GO_Env, inherits = FALSE)) {
            return(FALSE)
        }
    }
    isTRUE(org == get("hybrid_organism", envir = GO_Env)) &&
        isTRUE(keytype == get("hybrid_keytype", envir = GO_Env)) &&
        isTRUE(ont == get("hybrid_ont", envir = GO_Env))
}

#' Build a gene-to-GO annotation table restricted to given evidence codes
#'
#' Enrichment normally uses every annotation in the OrgDb, including
#' electronically inferred ones. Restricting the annotation by evidence code is
#' a methodological choice — for example dropping `IEP` when the gene sets came
#' from co-expression clustering, to avoid reasoning in a circle, or dropping
#' `IEA` to keep only manually reviewed annotations (#160).
#'
#' Note the consequence: genes whose annotations are all filtered away leave the
#' background as well, so the universe becomes "genes with at least one
#' annotation carrying an allowed evidence code". That is usually what is wanted,
#' but it does change the p-values beyond simply removing terms.
#'
#' @param OrgDb a loaded OrgDb object
#' @param keytype gene id type to report
#' @param evidence character vector of evidence codes to keep
#' @return a data.frame with the gene id, `GOALL` and `ONTOLOGYALL`
#' @importFrom AnnotationDbi keys
#' @noRd
go_annotation_by_evidence <- function(OrgDb, keytype, evidence) {
    evidence <- unique(as.character(evidence))
    if (!length(evidence)) {
        stop("`evidence` must contain at least one evidence code")
    }

    genes <- AnnotationDbi::keys(OrgDb, keytype = keytype)
    anno <- suppressMessages(
        AnnotationDbi::select(
            OrgDb,
            keys = genes,
            columns = c("GOALL", "EVIDENCEALL", "ONTOLOGYALL"),
            keytype = keytype
        )
    )

    anno <- anno[!is.na(anno$GOALL) & !is.na(anno$EVIDENCEALL), , drop = FALSE]
    anno <- anno[anno$EVIDENCEALL %in% evidence, , drop = FALSE]

    if (nrow(anno) == 0) {
        stop("no GO annotation left after keeping evidence code(s): ",
             paste(evidence, collapse = ", "),
             "; check `AnnotationDbi::columns(OrgDb)` for the codes available")
    }

    unique(anno[, c(keytype, "GOALL", "ONTOLOGYALL")])
}

#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi keytypes
#' @importFrom AnnotationDbi toTable
#' @importFrom GO.db GOTERM
get_GO_data <- function(OrgDb, ont, keytype, evidence = NULL) {
    GO_Env <- get_GO_Env()
    use_cached <- FALSE

    ## An evidence-filtered annotation is deliberately never cached: the cache is
    ## keyed by (organism, keytype) only, so storing a filtered annotation
    ## under that key — or serving the unfiltered one for a filtered request —
    ## would quietly answer with the wrong annotation.
    if (is.null(evidence) &&
        go_cache_usable(GO_Env, get_organism(OrgDb), ont, keytype)) {
        use_cached <- TRUE
    }

    if (use_cached) {
        goAnno <- get("goAnno", envir = GO_Env)
    } else {
        OrgDb <- load_OrgDb(OrgDb)
        kt <- keytypes(OrgDb)
        if (!keytype %in% kt) {
            stop("keytype is not supported...")
        }

        if (is.null(evidence)) {
        ## Build the all-ontology annotation: it is a superset of every
        ## single-ontology view, so the cache can serve any later `ont`
        ## without another mapIds pass. Filtering to the requested ontology
        ## happens below, after the cache lookup.
        goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
        go2gene <- suppressMessages(
            AnnotationDbi::mapIds(
                OrgDb,
                keys = names(goterms),
                column = keytype,
                keytype = "GOALL",
                multiVals = 'list'
            )
        )
        goAnno <- stack(go2gene)
        colnames(goAnno) <- c(keytype, "GOALL")
        goAnno <- unique(goAnno[!is.na(goAnno[, 1]), ])
        goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]

        assign("goAnno", goAnno, envir = GO_Env)
        assign("keytype", keytype, envir = GO_Env)
        assign("organism", get_organism(OrgDb), envir = GO_Env)
        } else {
            goAnno <- go_annotation_by_evidence(OrgDb, keytype, evidence)
        }
    }

    # Filter if needed (if cached was ALL but we want specific)
    if (ont != "ALL") {
        goAnno <- goAnno[goAnno$ONTOLOGYALL == ont, ]
    }

    # Build GSON
    GO2GENE <- unique(goAnno[, c(2, 1)])
    gsid2gene <- GO2GENE
    colnames(gsid2gene) <- c("gsid", "gene")
    
    termmap <- get_GO2TERM_table()
    gsid2name <- termmap[termmap$go_id %in% gsid2gene$gsid, ]
    colnames(gsid2name) <- c("gsid", "name")
    
    gson_obj <- gson::gson(gsid2gene = gsid2gene,
                           gsid2name = gsid2name,
                           species = get_organism(OrgDb),
                           gsname = paste0("GO_", ont),
                           keytype = keytype,
                           version = "unknown",
                           accessed_date = as.character(Sys.Date()))
    
    return(gson_obj)
}

get_GO_Env <- function() {
    if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos)
        assign(".GO_clusterProfiler_Env", new.env(), envir = envir)
    }
    get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}


## #' @importMethodsFrom AnnotationDbi Ontology
## #' @importFrom GO.db GOTERM
## #' @importMethodsFrom AnnotationDbi mappedkeys
## #' @importFrom plyr dlply
## #' @importFrom plyr .
## #' @importClassesFrom methods data.frame
## EXTID2TERMID.GO <- function(gene, ont, organism) {

##     gene <- as.character(gene)

##     ## get all goterms within the specific ontology
##     goterms <- Ontology(GOTERM)
##     goterms <- names(goterms[goterms == ont])

##     supported_Org <- getSupported_Org()
##     if (organism %in% supported_Org) {
##         mappedDb <- getGO2ALLEG_MappedDb(organism)

##         orgTerm <- mappedkeys(mappedDb)

##         ## narrow down goterms to specific organism
##         Terms <- goterms[goterms %in% orgTerm]

##         ## mapping GO to External gene ID
##         class(Terms) <- ont
##         GO2ExtID <- TERMID2EXTID(Terms, organism)

##         qGO2ExtID = lapply(GO2ExtID, function(i) gene[gene %in% i])
##         len <- sapply(qGO2ExtID, length)
##         notZero.idx <- len != 0
##         qGO2ExtID <- qGO2ExtID[notZero.idx]

##         len <- sapply(qGO2ExtID, length)
##         qGO2ExtID.df <- data.frame(GO=rep(names(qGO2ExtID), times=len),
##                                    ExtID=unlist(qGO2ExtID))

##         ExtID <- NULL ## to satisfy codetools
##         qExtID2GO <- dlply(qGO2ExtID.df, .(ExtID), function(i) as.character(i$GO))
##     } else {
##         oldwd <- getwd()

##         if (file.exists("EG2ALLGO.rda")) {
##             EG2ALLGO <- NULL # to satisfy codetools
##             load("EG2ALLGO.rda")
##             qExtID2GO <- EG2ALLGO[gene]
##             qExtID2GO <- lapply(qExtID2GO, function(i) i[i %in% goterms])
##         } else if (organism == "D39") {
##             dir <- system.file("extdata/D39/", package="clusterProfiler")
##             setwd(dir)
##         } else if (organism == "M5005") {
##             dir <- system.file("extdata/M5005/", package="clusterProfiler")
##             setwd(dir)
##         } else {
##             setwd(oldwd)
##             stop("GO mapping files not found in the working directory")
##         }

##         setwd(oldwd)
##     }
##     return(qExtID2GO)
## }

## #' @importMethodsFrom AnnotationDbi mget
## #' @importFrom GOSemSim getSupported_Org
## TERMID2EXTID.GO <- function(term, organism, ...) {
##     term <- as.character(term)

##     GO2ALLEG <- GO2EXTID(organism)
##     if (is(GO2ALLEG, "Go3AnnDbBimap")) {
##         GO2ExtID <- mget(term, GO2ALLEG, ifnotfound=NA)
##         GO2ExtID <- lapply(GO2ExtID, function(i) unique(i))
##     } else {
##         GO2ExtID <- GO2ALLEG[term]
##     }
##     return(GO2ExtID)
## }

## GO2EXTID <- function(organism) {
##     supported_Org <- getSupported_Org()
##     if (organism %in% supported_Org) {
##         GO2ALLEG <- getGO2ALLEG_MappedDb(organism)
##     } else {
##         oldwd <- getwd()
##         if(organism == "D39") {
##             dir <- system.file("extdata/D39/", package="clusterProfiler")
##             setwd(dir)
##         }
##         if(organism == "M5005") {
##             dir <- system.file("extdata/M5005/", package="clusterProfiler")
##             setwd(dir)
##         }
##         if (file.exists("GO2ALLEG.rda")) {
##             GO2ALLEG <- NULL # to satisfy codetools
##             load("GO2ALLEG.rda")
##         } else {
##             setwd(oldwd)
##             stop("GO Mapping file not found in the working directory")
##         }
##         setwd(oldwd)
##     }
##     return(GO2ALLEG)
## }

## #' @importMethodsFrom AnnotationDbi mappedkeys
## #' @importFrom GOSemSim getSupported_Org
## ALLEXTID.GO <- function(organism) {
##     supported_Org <- getSupported_Org()
##     if (organism %in% supported_Org) {
##         mappedDb <- getEG2GO_MappedDb(organism)
##         extID <- mappedkeys(mappedDb)
##     } else {
##         oldwd <- getwd()
##         if(organism == "D39") {
##             dir <- system.file("extdata/D39/", package="clusterProfiler")
##             setwd(dir)
##         }
##         if(organism == "M5005") {
##             dir <- system.file("extdata/M5005/", package="clusterProfiler")
##             setwd(dir)
##         }
##         if (file.exists("EG2ALLGO.rda")) {
##             EG2ALLGO <- NULL ## to satisfy codetools
##             load("EG2ALLGO.rda")
##             extID <- names(EG2ALLGO)
##         } else {
##             setwd(oldwd)
##             stop("GO mapping file not found in the working directory")
##         }
##         setwd(oldwd)
##     }
##     return(extID)
## }

#' drop GO term of specific level or specific terms (mostly too general).
#'
#'
#' @title dropGO
#' @param x an instance of 'enrichResult' or 'compareClusterResult'
#' @param level GO level
#' @param term GO term
#' @return modified version of x
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Ontology
#' @export
#' @author Guangchuang Yu
dropGO <- function(x, level = NULL, term = NULL) {
    if (!(is(x, "enrichResult") || is(x, "compareClusterResult"))) {
        stop(
            "x should be an instance of 'enrichResult' or 'compareClusterResult' ..."
        )
    }

    res <- as.data.frame(x)

    if (!is.null(level)) {
        if (is(x, "enrichResult")) {
            ont <- x@ontology
        } else {
            ont <- x@.call$ont
            if (is.null(ont)) {
                ## should be "MF", default value of enrichGO
                ## it's safe to determine from the output
                ont <- res$ID[1] %>% GOTERM[[.]] %>% Ontology
            }
        }

        tt <- getGOLevel(ont, level)
        term <- c(term, tt) %>% unique
    }
    if (is.null(term)) {
        return(x)
    }

    if (is(x, "enrichResult")) {
        res <- res[!res$ID %in% term, ]
        x@result <- res
    } else {
        res <- res[!res$ID %in% term, ]
        x@compareClusterResult <- res
    }

    return(x)
}
