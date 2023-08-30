getPPI_oldversion <- function(x, ID=1, taxID = "auto", limit = NULL, output = 'igraph') {
    output <- match.arg(output, c("igraph", "data.frame"))

    if (taxID == "auto") {
        if (is.null(x@organism) || length(x@organism) == 0) {
            stop("Unable to determine taxonomy ID. You can use the `getTaxID()` to get the taxonomy ID from scientific species name.\n")
        } 
        taxID <- getTaxID(x@organism)
    }

    genes <- paste(unlist(geneInCategory(x)[x$ID[ID]]), collapse = "%0d")

    if (is.null(limit)) {
        ## only use the genes in interest
        interaction_type <- "network?identifiers="
        limit <- ""
    } else {
        ## incorporate partner genes, size restricted by limit
        interaction_type <- "interaction_partners?identifiers="
        limit <- paste("&limit=",limit,sep = "")
    }

    url <- paste0("https://string-db.org/api/json/",
                interaction_type, genes,
                "&species=", taxID,
                limit
            )
   
    rlang::check_installed('jsonlite', 'for `getPPI_olderversion`.')

    res <- jsonlite::fromJSON(url)
    if (output == "data.frame") {
        return(res)
    }

    node <- unique(c(res$preferredName_A, res$preferredName_B))
    
    igraph::graph_from_data_frame(
        d = res[,c(3,4,6)],
        vertices = node,
        directed=F
    )    
}


#' networkParamsParser
#' 
#' parameters parser for [Getting the STRING network interactions](https://string-db.org/cgi/help.pl?sessionId=btsvnCeNrBk7).
#'
#' @param identifiers required parameter for multiple items, e.g. `c("PTCH1", "TP53", "BRCA1", "BRCA2")`
#' @param species NCBI taxon identifiers (e.g. Human is 9606, see: [STRING organisms](https://string-db.org/cgi/input.pl?input_page_active_form=organisms).
#' @param required_score threshold of significance to include a interaction, a number between 0 and 1000 (default depends on the network)
#' @param network_type network type: functional (default), physical
#' @param add_nodes adds a number of proteins with to the network based on their confidence score (default:1)
#' @param show_query_node_labels when available use submitted names in the preferredName column when (0 or 1) (default:0)
#' @param caller_identity your identifier for use.
#' 
#' @noRd
#' @return a list contain parameters for query
networkParamsParser <- function(
    identifiers,
    species,
    required_score = NULL,
    network_type = "functional",
    add_nodes = 1,
    show_query_node_labels = 0,
    caller_identity = NULL
) {
  # Format the identifiers
  identifiers <- paste(identifiers, collapse = "\n")
  
  # Check parameters
  if (missing(species)) {
    stop("Please provide an NCBI taxon identifier for the species.")
  }
  
  # Create parameters list
  params <- list(
    identifiers = identifiers,
    species = species,
    required_score = required_score,
    network_type = network_type,
    add_nodes = add_nodes,
    show_query_node_labels = show_query_node_labels,
    caller_identity = caller_identity
  )
  
  # Remove NULL elements from the list
  filtered_params <- Filter(Negate(is.null), params)
  return(filtered_params)
}

#' getPPI
#' 
#' [Getting the STRING network interactions](https://string-db.org/cgi/help.pl?sessionId=btsvnCeNrBk7).
#' 
#' @title getPPI
#' @param x an `enrichResult`` object or a vector of proteins, e.g. `c("PTCH1", "TP53", "BRCA1", "BRCA2")`
#' @param ID ID or index to extract genes in the enriched term(s) if `x` is an `enrichResult` object
#' @param taxID NCBI taxon identifiers (e.g. Human is 9606, see: [STRING organisms](https://string-db.org/cgi/input.pl?input_page_active_form=organisms).
#' @param required_score threshold of significance to include a interaction, a number between 0 and 1000 (default depends on the network)
#' @param network_type network type: functional (default), physical
#' @param add_nodes adds a number of proteins with to the network based on their confidence score (default:1)
#' @param show_query_node_labels when available use submitted names in the preferredName column when (0 or 1) (default:0)
#' @param output one of `data.frame` or `igraph`
#' @author Yonghe Xia and modified by Guangchuang Yu
#' @return a `data.frame` or an `igraph` object
#' @importFrom yulab.utils yread_tsv
#' @export
getPPI <- function(
    x, 
    ID=1, 
    taxID = "auto",     
    required_score = NULL,
    network_type = "functional",
    add_nodes = 0,
    show_query_node_labels = 0,
    output = 'igraph') {

    output <- match.arg(output, c("igraph", "data.frame"))
    if (taxID == "auto") {
        if (!inherits(x, 'enrichResult') || is.null(x@organism) || length(x@organism) == 0) {
            stop("Unable to determine taxonomy ID. You can use the `getTaxID()` to get the taxonomy ID from scientific species name.\n")
        } 
        taxID <- getTaxID(x@organism)
    }
    if (inherits(x, 'enrichResult')) {
        if (is(ID, 'numeric')) {
            id <- x$ID[ID]
        } else {
            id <- ID
        }
        genes <- unlist(geneInCategory(x)[id])
    } else {
        genes <- x
    }

    networkParams <- networkParamsParser(
        identifiers = genes,
        species = taxID,
        required_score = required_score,
        network_type = network_type,
        add_nodes = add_nodes,
        show_query_node_labels = show_query_node_labels
    )
    
    # Set stringDB base URL
    address <- "https://string-db.org"

    # Validate the address
    # 
    httr::stop_for_status(httr::GET(address))

    # read data from stringDB api
    response <- httr::GET(paste(address, "/api/tsv/network", sep = ""), query = networkParams)
    res <- yread_tsv(response$url, params = list(header = TRUE))

    if (output == "data.frame") {
        return(res)
    }

    node <- unique(c(res$preferredName_A, res$preferredName_B))
    
    igraph::graph_from_data_frame(
        d = unique(res[,c(3,4,6)]),
        vertices = node,
        directed=FALSE
    )    
}

