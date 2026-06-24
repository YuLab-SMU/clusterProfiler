#' Network-based Gene Set Enrichment Analysis of Gene Ontology
#'
#' @title nseGO
#' @param geneList order ranked geneList
#' @param network edge list (data.frame with 2 or 3 columns) or sparse matrix
#' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param OrgDb OrgDb
#' @param keyType keytype of gene
#' @param mode propagation mode, one of "evidence" and "signed"
#' @param p restart probability for RWR
#' @param specific_weight logical, whether to apply gene specificity weighting
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param threshold convergence threshold for RWR
#' @param maxIter maximal number of RWR iterations
#' @param verbose print message or not
#' @param ... other parameters passed to \code{enrichit::nsea_gson()}
#' @importClassesFrom enrichit nseaResult
#' @return nseaResult object
#' @export
nseGO <- function(geneList,
                  network,
                  ont = "BP",
                  OrgDb,
                  keyType = "ENTREZID",
                  mode = c("evidence", "signed"),
                  p = 0.5,
                  specific_weight = FALSE,
                  minGSSize = 10,
                  maxGSSize = 500,
                  threshold = 1e-9,
                  maxIter = 100,
                  verbose = TRUE,
                  ...) {
    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

    GO_DATA <- get_GO_data(OrgDb, ont, keyType)

    res <- enrichit::nsea_gson(
        geneList = geneList,
        network = network,
        gson = GO_DATA,
        mode = mode,
        p = p,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_go_result(res, GO_DATA, ont, OrgDb, keyType)
}

#' Multi-layer Network-based Gene Set Enrichment Analysis of Gene Ontology
#'
#' @title mnseGO
#' @param seed_list named list of named numeric vectors, one per layer
#' @param networks named list of layer-specific networks
#' @param couplings data.frame of inter-layer edges
#' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param OrgDb OrgDb
#' @param keyType keytype of gene
#' @param mode propagation mode, one of "evidence" and "signed"
#' @param layer_weights optional named numeric vector of layer weights
#' @param collapse one of "weighted_mean", "sum", "mean", or "max_abs"
#' @param target_layer optional layer name to export scores from
#' @param output_space one of "union" and "gene"
#' @param p restart probability for RWR
#' @param interlayer_strength global scaling factor for coupling edges
#' @param specific_weight logical, whether to apply gene specificity weighting
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param threshold convergence threshold for RWR
#' @param maxIter maximal number of RWR iterations
#' @param verbose print message or not
#' @param ... other parameters passed to \code{enrichit::mnsea_gson()}
#' @importClassesFrom enrichit mnseaResult
#' @return mnseaResult object
#' @export
mnseGO <- function(seed_list,
                   networks,
                   couplings,
                   ont = "BP",
                   OrgDb,
                   keyType = "ENTREZID",
                   mode = c("evidence", "signed"),
                   layer_weights = NULL,
                   collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                   target_layer = NULL,
                   output_space = c("union", "gene"),
                   p = 0.5,
                   interlayer_strength = 1,
                   specific_weight = FALSE,
                   minGSSize = 10,
                   maxGSSize = 500,
                   threshold = 1e-9,
                   maxIter = 100,
                   verbose = TRUE,
                   ...) {
    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)

    GO_DATA <- get_GO_data(OrgDb, ont, keyType)

    res <- enrichit::mnsea_gson(
        seed_list = seed_list,
        networks = networks,
        couplings = couplings,
        gson = GO_DATA,
        mode = mode,
        layer_weights = layer_weights,
        collapse = collapse,
        target_layer = target_layer,
        output_space = output_space,
        p = p,
        interlayer_strength = interlayer_strength,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_go_result(res, GO_DATA, ont, OrgDb, keyType)
}

#' Network-based Gene Set Enrichment Analysis of KEGG
#'
#' @title nseKEGG
#' @param geneList order ranked geneList
#' @param network edge list (data.frame with 2 or 3 columns) or sparse matrix
#' @param organism supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html',
#'   or a GSON object
#' @param keyType one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
#' @param mode propagation mode, one of "evidence" and "signed"
#' @param p restart probability for RWR
#' @param specific_weight logical, whether to apply gene specificity weighting
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param threshold convergence threshold for RWR
#' @param maxIter maximal number of RWR iterations
#' @param verbose print message or not
#' @param use_internal_data logical, use KEGG.db or latest online KEGG data
#' @param ... other parameters passed to \code{enrichit::nsea_gson()}
#' @importClassesFrom enrichit nseaResult
#' @return nseaResult object
#' @export
nseKEGG <- function(geneList,
                    network,
                    organism = "hsa",
                    keyType = "kegg",
                    mode = c("evidence", "signed"),
                    p = 0.5,
                    specific_weight = FALSE,
                    minGSSize = 10,
                    maxGSSize = 500,
                    threshold = 1e-9,
                    maxIter = 100,
                    verbose = TRUE,
                    use_internal_data = FALSE,
                    ...) {
    kegg_info <- .prepare_nse_kegg_data(organism, keyType, use_internal_data)

    res <- enrichit::nsea_gson(
        geneList = geneList,
        network = network,
        gson = kegg_info$gson,
        mode = mode,
        p = p,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_kegg_result(res, kegg_info$species, kegg_info$keyType)
}

#' Network-based Gene Set Enrichment Analysis of KEGG Module
#'
#' @title nseMKEGG
#' @inheritParams nseKEGG
#' @param ... other parameters passed to \code{enrichit::nsea_gson()}
#' @return nseaResult object
#' @export
nseMKEGG <- function(geneList,
                     network,
                     organism = "hsa",
                     keyType = "kegg",
                     mode = c("evidence", "signed"),
                     p = 0.5,
                     specific_weight = FALSE,
                     minGSSize = 10,
                     maxGSSize = 500,
                     threshold = 1e-9,
                     maxIter = 100,
                     verbose = TRUE,
                     ...) {
    kegg_info <- .prepare_nse_kegg_data(organism, keyType, FALSE, "MKEGG")

    res <- enrichit::nsea_gson(
        geneList = geneList,
        network = network,
        gson = kegg_info$gson,
        mode = mode,
        p = p,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_kegg_result(
        res,
        species = kegg_info$species,
        keyType = kegg_info$keyType,
        setType = "MKEGG",
        append_category = FALSE
    )
}

#' Multi-layer Network-based Gene Set Enrichment Analysis of KEGG
#'
#' @title mnseKEGG
#' @param seed_list named list of named numeric vectors, one per layer
#' @param networks named list of layer-specific networks
#' @param couplings data.frame of inter-layer edges
#' @param organism supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html',
#'   or a GSON object
#' @param keyType one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
#' @param mode propagation mode, one of "evidence" and "signed"
#' @param layer_weights optional named numeric vector of layer weights
#' @param collapse one of "weighted_mean", "sum", "mean", or "max_abs"
#' @param target_layer optional layer name to export scores from
#' @param output_space one of "union" and "gene"
#' @param p restart probability for RWR
#' @param interlayer_strength global scaling factor for coupling edges
#' @param specific_weight logical, whether to apply gene specificity weighting
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param threshold convergence threshold for RWR
#' @param maxIter maximal number of RWR iterations
#' @param verbose print message or not
#' @param use_internal_data logical, use KEGG.db or latest online KEGG data
#' @param ... other parameters passed to \code{enrichit::mnsea_gson()}
#' @importClassesFrom enrichit mnseaResult
#' @return mnseaResult object
#' @export
mnseKEGG <- function(seed_list,
                     networks,
                     couplings,
                     organism = "hsa",
                     keyType = "kegg",
                     mode = c("evidence", "signed"),
                     layer_weights = NULL,
                     collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                     target_layer = NULL,
                     output_space = c("union", "gene"),
                     p = 0.5,
                     interlayer_strength = 1,
                     specific_weight = FALSE,
                     minGSSize = 10,
                     maxGSSize = 500,
                     threshold = 1e-9,
                     maxIter = 100,
                     verbose = TRUE,
                     use_internal_data = FALSE,
                     ...) {
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)
    kegg_info <- .prepare_nse_kegg_data(organism, keyType, use_internal_data)

    res <- enrichit::mnsea_gson(
        seed_list = seed_list,
        networks = networks,
        couplings = couplings,
        gson = kegg_info$gson,
        mode = mode,
        layer_weights = layer_weights,
        collapse = collapse,
        target_layer = target_layer,
        output_space = output_space,
        p = p,
        interlayer_strength = interlayer_strength,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_kegg_result(res, kegg_info$species, kegg_info$keyType)
}

#' Multi-layer Network-based Gene Set Enrichment Analysis of KEGG Module
#'
#' @title mnseMKEGG
#' @inheritParams mnseKEGG
#' @param ... other parameters passed to \code{enrichit::mnsea_gson()}
#' @return mnseaResult object
#' @export
mnseMKEGG <- function(seed_list,
                      networks,
                      couplings,
                      organism = "hsa",
                      keyType = "kegg",
                      mode = c("evidence", "signed"),
                      layer_weights = NULL,
                      collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                      target_layer = NULL,
                      output_space = c("union", "gene"),
                      p = 0.5,
                      interlayer_strength = 1,
                      specific_weight = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500,
                      threshold = 1e-9,
                      maxIter = 100,
                      verbose = TRUE,
                      ...) {
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)
    kegg_info <- .prepare_nse_kegg_data(organism, keyType, FALSE, "MKEGG")

    res <- enrichit::mnsea_gson(
        seed_list = seed_list,
        networks = networks,
        couplings = couplings,
        gson = kegg_info$gson,
        mode = mode,
        layer_weights = layer_weights,
        collapse = collapse,
        target_layer = target_layer,
        output_space = output_space,
        p = p,
        interlayer_strength = interlayer_strength,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_kegg_result(
        res,
        species = kegg_info$species,
        keyType = kegg_info$keyType,
        setType = "MKEGG",
        append_category = FALSE
    )
}

#' Network-based Gene Set Enrichment Analysis of WikiPathways
#'
#' @title nseWP
#' @param geneList order ranked geneList
#' @param network edge list (data.frame with 2 or 3 columns) or sparse matrix
#' @param organism supported organisms listed by \code{get_wp_organisms()}, or a GSON object
#' @param mode propagation mode, one of "evidence" and "signed"
#' @param p restart probability for RWR
#' @param specific_weight logical, whether to apply gene specificity weighting
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param threshold convergence threshold for RWR
#' @param maxIter maximal number of RWR iterations
#' @param verbose print message or not
#' @param ... other parameters passed to \code{enrichit::nsea_gson()}
#' @return nseaResult object
#' @export
nseWP <- function(geneList,
                  network,
                  organism,
                  mode = c("evidence", "signed"),
                  p = 0.5,
                  specific_weight = FALSE,
                  minGSSize = 10,
                  maxGSSize = 500,
                  threshold = 1e-9,
                  maxIter = 100,
                  verbose = TRUE,
                  ...) {
    wp_info <- .prepare_nse_wp_data(organism)

    res <- enrichit::nsea_gson(
        geneList = geneList,
        network = network,
        gson = wp_info$gson,
        mode = mode,
        p = p,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_wp_result(res, wp_info$species, wp_info$keyType)
}

#' Multi-layer Network-based Gene Set Enrichment Analysis of WikiPathways
#'
#' @title mnseWP
#' @param seed_list named list of named numeric vectors, one per layer
#' @param networks named list of layer-specific networks
#' @param couplings data.frame of inter-layer edges
#' @param organism supported organisms listed by \code{get_wp_organisms()}, or a GSON object
#' @param mode propagation mode, one of "evidence" and "signed"
#' @param layer_weights optional named numeric vector of layer weights
#' @param collapse one of "weighted_mean", "sum", "mean", or "max_abs"
#' @param target_layer optional layer name to export scores from
#' @param output_space one of "union" and "gene"
#' @param p restart probability for RWR
#' @param interlayer_strength global scaling factor for coupling edges
#' @param specific_weight logical, whether to apply gene specificity weighting
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param threshold convergence threshold for RWR
#' @param maxIter maximal number of RWR iterations
#' @param verbose print message or not
#' @param ... other parameters passed to \code{enrichit::mnsea_gson()}
#' @return mnseaResult object
#' @export
mnseWP <- function(seed_list,
                   networks,
                   couplings,
                   organism,
                   mode = c("evidence", "signed"),
                   layer_weights = NULL,
                   collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                   target_layer = NULL,
                   output_space = c("union", "gene"),
                   p = 0.5,
                   interlayer_strength = 1,
                   specific_weight = FALSE,
                   minGSSize = 10,
                   maxGSSize = 500,
                   threshold = 1e-9,
                   maxIter = 100,
                   verbose = TRUE,
                   ...) {
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)
    wp_info <- .prepare_nse_wp_data(organism)

    res <- enrichit::mnsea_gson(
        seed_list = seed_list,
        networks = networks,
        couplings = couplings,
        gson = wp_info$gson,
        mode = mode,
        layer_weights = layer_weights,
        collapse = collapse,
        target_layer = target_layer,
        output_space = output_space,
        p = p,
        interlayer_strength = interlayer_strength,
        specific_weight = specific_weight,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        threshold = threshold,
        maxIter = maxIter,
        verbose = verbose,
        ...
    )

    .finalize_nse_wp_result(res, wp_info$species, wp_info$keyType)
}

.finalize_nse_go_result <- function(res, GO_DATA, ont, OrgDb, keyType) {
    if (is.null(res)) {
        return(res)
    }

    if (keyType == "SYMBOL") {
        res@readable <- TRUE
    }
    res@organism <- get_organism(OrgDb)
    res@setType <- ont
    res@keytype <- keyType

    if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA)
    }

    res
}

.prepare_nse_kegg_data <- function(organism,
                                   keyType,
                                   use_internal_data = FALSE,
                                   kegg_type = "KEGG") {
    if (inherits(organism, "character") && organism == "cpd") {
        organism <- gson_cpd()
    }

    if (inherits(organism, "character")) {
        species <- organismMapper(organism)
        if (use_internal_data) {
            KEGG_DATA <- get_data_from_KEGG_db(species)
        } else {
            KEGG_DATA <- prepare_KEGG(species, kegg_type, keyType)
        }
    } else if (inherits(organism, "GSON")) {
        KEGG_DATA <- organism
        species <- KEGG_DATA@species
        keyType <- KEGG_DATA@keytype
    } else {
        stop("organism should be a species name or a GSON object")
    }

    list(gson = KEGG_DATA, species = species, keyType = keyType)
}

.prepare_nse_wp_data <- function(organism) {
    if (inherits(organism, "character")) {
        wp_data <- gson_WP(organism)
        species <- organism
        keyType <- "ENTREZID"
    } else if (inherits(organism, "GSON")) {
        wp_data <- organism
        species <- wp_data@species
        keyType <- wp_data@keytype
        if (is.null(keyType) || identical(keyType, "")) {
            keyType <- "ENTREZID"
        }
    } else {
        stop("organism should be a species name or a GSON object")
    }

    list(gson = wp_data, species = species, keyType = keyType)
}

.finalize_nse_kegg_result <- function(res,
                                      species,
                                      keyType,
                                      setType = "KEGG",
                                      append_category = TRUE) {
    if (is.null(res)) {
        return(res)
    }

    res@organism <- species
    res@setType <- setType
    res@keytype <- keyType
    if (append_category) {
        res <- append_kegg_category(res)
    }
    res
}

.finalize_nse_wp_result <- function(res, species, keyType) {
    if (is.null(res)) {
        return(res)
    }

    res@organism <- species
    res@setType <- "WikiPathways"
    res@keytype <- keyType
    res
}
