##' download the latest version of KEGG pathway and stored in a 'GSON' object
##'
##'
##' @title gson_KEGG
##' @param species species
##' @param KEGG_Type one of "KEGG" and "MKEGG"
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'.
##' @return a 'GSON' object
##' @author Guangchuang Yu
##' @importFrom gson gson
##' @export
gson_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    x <- download_KEGG(species, KEGG_Type, keyType)
    gsid2gene <- setNames(x[[1]], c("gsid", "gene"))
    gsid2name <- setNames(x[[2]], c("gsid", "name"))
    version <- kegg_release(species)
    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name,
        species = species,
        gsname = "KEGG",
        version = version,
        accessed_date = as.character(Sys.Date(),
        keytype = "ENTREZID")
    )
}

##' download the latest version of KEGG pathway and stored in a 'GSON' object
##'
##'
##' @title gson_KEGG
##' @param OrgDb OrgDb
##' @param keytype keytype of genes.
##' @param ont one of "BP", "MF", "CC", and "ALL"
##' @return a 'GSON' object
##' @importFrom gson gson
##' @export
gson_GO <- function(OrgDb, keytype = 'ENTREZID', ont = "BP") {

    if (is(OrgDb, "character")) {
        require(OrgDb, character.only = TRUE)
        OrgDb <- eval(parse(text = OrgDb))
    }

    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    if (ont != "ALL") {
        goterms <- goterms[goterms == ont]
    }
    go2gene <- suppressMessages(
        AnnotationDbi::mapIds(OrgDb, keys=names(goterms), column=keytype,
                            keytype="GOALL", multiVals='list')
    )
    goAnno <- stack(go2gene)   
    gsid2gene <- goAnno[, c(2,1)]
    colnames(gsid2gene) <- c("gsid", "gene")
    gsid2gene <- unique(gsid2gene[!is.na(gsid2gene[,2]), ]) 

    termname <- AnnotationDbi::Term(GO.db::GOTERM)
    gsid2name <- data.frame(gsid = names(termname),
                            name = termname)
    species <- AnnotationDbi::species(OrgDb)
    m <- AnnotationDbi::metadata(OrgDb)
    version <- m$value[m$name == "GOSOURCEDATE"]
    # gsname <- m$value[m$name == 'GOSOURCENAME']
    gsname <- paste(m$value[m$name == 'GOSOURCENAME'], ont, sep = ";")

   # gene2name
    genes <- unique(gsid2gene[, 2])
    gene2name <- bitr(geneID = genes, fromType = keytype, 
        toType = "SYMBOL", OrgDb = OrgDb, drop = TRUE)


    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name,
        gene2name = gene2name,
        species = species,
        gsname = gsname,
        version = version,
        accessed_date = as.character(Sys.Date()),
        keytype = keytype
    )
}
##' Download the latest version of WikiPathways data and stored in a 'GSON' object
##'
##'
##' @title gson_WP
##' @param organism supported organism, which can be accessed via the get_wp_organisms() function.
##' @export
gson_WP <- function(organism) {
    get_wp_data(organism, output = "GSON")
}


#' Build a gson object that annotate Gene Ontology
#'
#' @param data a two-column data.frame of original GO annotation. The columns are "gene_id" and "go_id".
#' @param ont type of GO annotation, which is "ALL", "BP", "MF", or "CC". default: "ALL".
#' @param species name of species. Default: NULL.
#' @param ... pass to `gson::gson()` constructor.
#'
#' @return a `gson` instance
#' @export
#'
#' @examples
#'  data = data.frame(gene_id = "gene1", 
#'                    go_id = c("GO:0035492", "GO:0009764", "GO:0031063", "GO:0033714", "GO:0036349"))
#'  gson_go_mapper(data, species = "E. coli")
gson_GO_mapper = function(data, 
                     ont = c("ALL", "BP", "CC", "MF"), 
                     species = NULL, 
                     ...){
  ont = match.arg(ont)
  
  data = unique(data) # cleanup
  if (nrow(data) == 0) {
    simpleError("Data is empty in this call.")
  }
  
  # resources from `GO.db`
  goterms = AnnotationDbi::Ontology(GO.db::GOTERM)
  termname = AnnotationDbi::Term(GO.db::GOTERM)
  go.db_info = GO.db::GO_dbInfo()
  go.db_source_date = go.db_info[go.db_info$name == "GOSOURCEDATE", "value"]
  ancestor_map = lapply(c(GO.db::GOBPANCESTOR, 
                          GO.db::GOCCANCESTOR,
                          GO.db::GOMFANCESTOR),
                        as.list) %>%
    unlist(recursive = FALSE)
  
  # filter GO terms
  data[["ontology"]] = goterms[data[["go_id"]]]
  n_na_ont = sum(is.na(data[["ontology"]]))
  if ( n_na_ont > 0){
    warning(sprintf("%s GO term(s) are too new for current `GO.db` [source date: %s],\n  and are to be dropped. Consider to update `GO.db` if possible.",
                    n_na_ont, 
                    go.db_source_date))
  }
  
  # map to GO ancestor
  ancestor_list = ancestor_map[data$go_id]
  names(ancestor_list) = data$gene_id
  ancestor_go = AnnotationDbi::unlist2(ancestor_list)
  
  # gsid2gene
  gsid2gene = data.frame(
    gsid = c(ancestor_go, data$go_id),
    gene = c(names(ancestor_go), data$gene_id),
    ontology = goterms[c(ancestor_go, data$go_id)]
  ) %>%
    dplyr::filter(.data$gsid != "all") %>%
    unique()
  
  if (ont != "ALL"){
    gsid2gene = gsid2gene %>%
      dplyr::filter(.data$ontology == ont)
  }
  
  # gsid2name
  uniq_gsid = unique(gsid2gene$gsid) %>% as.character()
  gsid2name = data.frame(
    gsid = uniq_gsid,
    name = termname[uniq_gsid] %>% as.character()
  )
  
  # construct `gson` object
  gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = species,
    gsname = paste0("Gene Ontology: ", ont),
    version = sprintf("[GO.db source date: %s]", go.db_source_date),
    accessed_date = as.character(Sys.Date()),
    ...
  )
}


#' Build KEGG annotation for novel species using KEGG Mapper
#' 
#' KEGG Mapper service can annotate protein sequences for novel species with KO database, 
#' and KO annotation need to be converted into Pathway or Module annotation, which
#' can then be used in `clusterProfiler`
#' 
#' @export
#' @return  a gson instance
#' @param file the name of the file which comes from the KEGG Mapper service, see Details for file format
#' @param format string indicate format of KEGG Mapper result
#' @param type string indicate annotation database
#' @param species your species, NULL if ignored
#' @param ... pass to gson::gson()
#' @details File is a two-column dataset with K numbers in the second column, optionally preceded by 
#'          the user's identifiers in the first column. This is consistent with the output 
#'          files of automatic annotation servers, BlastKOALA, GhostKOALA, and KofamKOALA. 
#'          KOALA (KEGG Orthology And Links Annotation) is KEGG's internal annotation tool 
#'          for K number assignment of KEGG GENES using SSEARCH computation. BlastKOALA 
#'          and GhostKOALA assign K numbers to the user's sequence data by BLAST and GHOSTX 
#'          searches, respectively, against a nonredundant set of KEGG GENES. KofamKOALA 
#'          is a new member of the KOALA family available at GenomeNet using the HMM profile 
#'          search, rather than the sequence similarity search, for K number assignment. 
#'          see https://www.kegg.jp/blastkoala/, https://www.kegg.jp/ghostkoala/ and 
#'          https://www.genome.jp/tools/kofamkoala/ for more information.
#' @examples 
#' \dontrun{
#'  file = system.file('extdata', "kegg_mapper_blast.txt", package='clusterProfiler')
#'  gson_KEGG_mapper(file, format = "BLAST", type = "pathway")
#' }
gson_KEGG_mapper = function(file, 
                       format = c("BLAST","Ghost","Kofam"), 
                       type = c("pathway","module"),
                       species = NULL,
                       ...){
  message("Please be aware that this process need a active internet connection and may take a while to finish.")
  format = match.arg(format)
  type = match.arg(type)
  if (format %in% c("BLAST", "Ghost")){
    protein2ko = utils::read.delim(file, header = FALSE, col.names = c("id","ko")) %>% 
      dplyr::filter(.data[["ko"]] != "") %>%
      dplyr::mutate_all("as.character")
  } else if (format == "Kofam"){
    protein2ko = utils::read.delim(file, comment.char = "#", header = FALSE) %>% 
      tidyr::separate(1, into = c("star", 'id', 'ko'), sep = " +", extra = "drop") %>% 
      dplyr::select(-1) %>%
      tidyr::drop_na() %>%
      dplyr::mutate_all("as.character")
  } else {
    simpleError("Please specify a valid file format for your KEGG Mapper result.")
  }
  
  if (type == "pathway"){
    all_pathway = kegg_list("pathway")
    colnames(all_pathway) = c("pathway","pathway_name")
    all_pathway %<>%
      dplyr::mutate_at(.vars = "pathway", .funs = "remove_db_prefix")
    
    ko2pathway = kegg_link("pathway", "ko")
    colnames(ko2pathway) = c("ko", "pathway")
    ko2pathway %<>%
      dplyr::mutate_all(.funs = "remove_db_prefix") %>%
      # remove redundant ko prefix pathway ids, only keep those like map00010
      # dplyr::filter(stringr::str_starts(.data[["pathway"]], "map"))
      dplyr::filter(grepl(.data$pathway, pattern = "^map"))
      
    # some KOs don't have pathway mapping, in which case NA will be induced at pathway
    gsid2gene = protein2ko %>%
      dplyr::left_join(ko2pathway, by = "ko") %>%
      dplyr::select(c("pathway","id"))
    na_id = gsid2gene %>% 
      dplyr::filter(is.na(.data[["pathway"]])) %>%
      dplyr::pull(.data[["id"]])
    if (length(na_id) > 0){
      warning(sprintf("%d lines of KEGG mapper result didn't have pathway map, and they are ignored.\n  Totally affected %d unique gene ids.\n  The first five id are:\n%s.",
                      length(na_id), 
                      dplyr::n_distinct(na_id), 
                      paste0("    ", na_id[1:5], collapse = "\n")))
      gsid2gene = gsid2gene %>%
        dplyr::filter(!is.na(.data[["pathway"]]))
    }
    
    gsid2name = gsid2gene %>%
      dplyr::distinct(.data[["pathway"]]) %>%
      dplyr::left_join(all_pathway, by = "pathway")
  } else if (type == "module"){
    # all kegg module
    all_module = kegg_list("module")
    colnames(all_module) = c("module","module_name")
    all_module = all_module %>%
      dplyr::mutate_at(.vars = "module", .funs = "remove_db_prefix")
    
    # ko2module
    ko2module = kegg_link("module", "ko")
    colnames(ko2module) = c( "ko", "module")
    ko2module = ko2module %>%
      dplyr::mutate_all(.funs = "remove_db_prefix")
    
    # module to gene
    gsid2gene = protein2ko %>%
      dplyr::left_join(ko2module, by = "ko") %>%
      dplyr::select(c("module", "id"))
    
    # similar to pathway, some KOs don't mapp to any module
    na_id = gsid2gene %>% 
      dplyr::filter(is.na(.data[["module"]])) %>%
      dplyr::pull(.data[["id"]])
    if (length(na_id) > 0){
      warning(sprintf("%d lines of KEGG mapper result didn't belong to any module, and they are ignored.\n  Totally affected %d unique gene ids.\n  The first five id are:\n%s.",
                      length(na_id), 
                      dplyr::n_distinct(na_id), 
                      paste0("    ", na_id[1:5], collapse = "\n")))
      gsid2gene = gsid2gene %>%
        dplyr::filter(!is.na(.data[["module"]]))
      
    }
    
    # gsid2name
    gsid2name = gsid2gene %>%
      dplyr::distinct(.data[["module"]]) %>%
      dplyr::left_join(all_module, by = "module")
  } else {
    simpleError('Please specify your target database (currently "pathway" or "module").')
  }
  
  # construct gson object
  release = kegg_release(type)
  colnames(gsid2gene) = c("gsid", "gene")
  colnames(gsid2name) = c("gsid", "name")
  gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = species,
    gsname = paste0("KEGG_", type),
    version = release,
    accessed_date = as.character(Sys.Date()),
    ...
  )
}

remove_db_prefix = function(x){
  gsub("^[a-z]+:", "", x)
}

kegg_release = function(db){
  url = paste0("https://rest.kegg.jp/info/", db)
  y = readLines(url)
  release <- sub("\\w+\\s+", "", y[grep('Release', y)])
  return(release)
}