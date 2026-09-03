#' @title Interpret Enrichment Results Using LLMs
#' @description
#' Functions for interpreting functional enrichment analysis results using
#' Large Language Models. Supports single-call interpretation, multi-agent
#' deep analysis, and hierarchical cluster strategies.
#'
#' Built on top of aisdk's `generate_object()` for reliable structured output,
#' and the Agent/Session system for multi-agent pipelines.
#' @name interpret
NULL

# ============================================================================
# Schema Definitions
# ============================================================================
#' @importFrom aisdk generate_object z_object z_string z_array z_enum ChatSession Agent
.interpretation_schema <- function() {
  z_object(
    overview = z_string("High-level summary of the key biological processes identified"),
    regulatory_drivers = z_string(
      "Key transcription factors or master regulators and their potential role in driving the observed pathways"
    ),
    key_mechanisms = z_string(
      "Explanation of underlying biological mechanisms, grouping pathways into major themes"
    ),
    crosstalk = z_string(
      "Discussion of potential interactions and regulatory networks between pathways"
    ),
    hypothesis = z_string(
      "A coherent biological hypothesis connecting the pathways (What) to biological meaning (So What)"
    ),
    narrative = z_string(
      "A cohesive paragraph suitable for the Results or Discussion section of a scientific paper"
    ),
    refined_network = z_array(
      z_object(
        source = z_string("Source gene or protein"),
        target = z_string("Target gene or protein"),
        interaction = z_string("Type of interaction: binding, activation, inhibition, etc."),
        reason = z_string("Why this edge is biologically relevant")
      ),
      description = "Core regulatory network edges"
    ),
    network_evidence = z_string(
      "Description of specific protein complexes or signaling modules that support the conclusion"
    )
  )
}

.annotation_schema <- function() {
  z_object(
    cell_type = z_string("The identified cell type label"),
    confidence = z_enum(c("High", "Medium", "Low"), description = "Confidence level of the annotation"),
    reasoning = z_string(
      "Explanation of why this cell type was assigned, citing specific markers or pathways from the input"
    ),
    regulatory_drivers = z_string("Key TFs or master regulators that define this cell type/state"),
    markers = z_array(
      z_string(),
      description = "Key markers or pathways from the input that support this decision"
    ),
    refined_network = z_array(
      z_object(
        source = z_string("Source gene or protein"),
        target = z_string("Target gene or protein"),
        interaction = z_string("Type of interaction"),
        reason = z_string("Why this edge is relevant")
      ),
      description = "Core regulatory network edges"
    ),
    network_evidence = z_string("Protein complexes or signaling modules that support the conclusion")
  )
}

.annotation_refinement_schema <- function() {
  z_object(
    cell_type = z_string("The final identified cell type label (refined or corrected)"),
    refinement_status = z_enum(
      c("Confirmed", "Refined", "Corrected"),
      description = "Whether the preliminary annotation was confirmed, refined, or corrected"
    ),
    confidence = z_enum(c("High", "Medium", "Low"), description = "Confidence level"),
    reasoning = z_string(
      "Why you confirmed, refined, or corrected the label, citing specific evidence"
    ),
    regulatory_drivers = z_string("Key TFs or master regulators"),
    markers = z_array(z_string(), description = "Supporting markers or pathways"),
    refined_network = z_array(
      z_object(
        source = z_string("Source gene or protein"),
        target = z_string("Target gene or protein"),
        interaction = z_string("Type of interaction"),
        reason = z_string("Why this edge is relevant")
      ),
      description = "Core regulatory network edges"
    ),
    network_evidence = z_string("Protein complexes or signaling modules supporting the conclusion")
  )
}

.phenotype_schema <- function() {
  z_object(
    phenotype = z_string("A concise label for the biological phenotype or state"),
    confidence = z_enum(c("High", "Medium", "Low"), description = "Confidence level"),
    reasoning = z_string("How the enriched terms support this phenotype"),
    regulatory_drivers = z_string("Key TFs or master regulators that drive this phenotype"),
    key_processes = z_array(
      z_string(),
      description = "Key pathways or terms that define this phenotype"
    ),
    refined_network = z_array(
      z_object(
        source = z_string("Source gene or protein"),
        target = z_string("Target gene or protein"),
        interaction = z_string("Type of interaction"),
        reason = z_string("Why this edge is relevant")
      ),
      description = "Core regulatory network edges"
    ),
    network_evidence = z_string("Protein complexes or signaling modules supporting the conclusion")
  )
}

.cleaner_schema <- function() {
  z_object(
    kept_pathways = z_array(z_string(), description = "Pathway names to keep"),
    discarded_pathways = z_array(z_string(), description = "Discarded pathway names"),
    reasoning = z_string("Explanation of the filtering strategy used")
  )
}

.detective_schema <- function() {
  z_object(
    key_drivers = z_array(
      z_string(),
      description = "Top 3-5 driver genes (TFs, Kinases)"
    ),
    functional_modules = z_array(
      z_string(),
      description = "Identified functional modules (e.g., TCR Complex, Cell Cycle G1/S)"
    ),
    refined_network = z_array(
      z_object(
        source = z_string("Source gene"),
        target = z_string("Target gene"),
        interaction = z_string("Type: activation, inhibition, binding, etc."),
        reason = z_string("Evidence for this interaction")
      ),
      description = "Core regulatory sub-network"
    ),
    network_evidence = z_string("How the network supports the identified drivers")
  )
}

# ============================================================================
# Model ID Helper
# ============================================================================

#' @title Infer Model ID
#' @description
#' Maps bare model names to the aisdk `provider:model` format for backward
#' compatibility. Emits a warning when guessing and suggests the explicit form.
#' If the model already contains a colon, it is returned as-is.
#'
#' @param model A model string, either bare (e.g., "deepseek-chat") or fully
#'   qualified (e.g., "deepseek:deepseek-chat").
#' @return A string in `provider:model` format.
#' @keywords internal
infer_model_id <- function(model) {
  if (!is.character(model) || length(model) != 1 || !nzchar(model)) {
    rlang::abort("`model` must be a non-empty model ID string.")
  }

  if (grepl(":", model, fixed = TRUE)) {
    return(model)
  }

  provider <- NULL

  patterns <- list(
    deepseek = "^deepseek-",
    openai = "^(gpt-|o[0-9]|chatgpt-)",
    anthropic = "claude",
    gemini = "^gemini-",
    bailian = "^(qwen-|qwq-|glm-)",
    volcengine = "^(doubao-|deepseek-r1)",
    stepfun = "^step-",
    xai = "^grok-",
    nvidia = "^(nvidia/|meta/llama)",
    openrouter = "/"
  )

  for (name in names(patterns)) {
    if (grepl(patterns[[name]], model, ignore.case = TRUE)) {
      provider <- name
      break
    }
  }

  if (is.null(provider)) {
    rlang::abort(c(
      paste0("Cannot infer provider for model: ", model),
      "i" = "Use the explicit 'provider:model' format.",
      "i" = "Example: 'gemini:gemini-3-flash-preview', 'openai:gpt-4o', 'deepseek:deepseek-chat'"
    ))
  }

  full_id <- paste0(provider, ":", model)
  rlang::warn(c(
    paste0("Inferred model ID as '", full_id, "'. Consider using the explicit format."),
    "i" = paste0("Use: model = '", full_id, "'")
  ))
  full_id
}

#' @keywords internal
.normalize_interpret_model <- function(model = NULL) {
  if (is.null(model) || inherits(model, "LanguageModelV1")) {
    return(model)
  }

  infer_model_id(model)
}

# ============================================================================
# Input Processing
# ============================================================================

#' @keywords internal
process_enrichment_input <- function(x, n_pathways) {
  get_df <- function(obj) {
    if (inherits(obj, "compareClusterResult") || inherits(obj, "enrichResult") || inherits(obj, "gseaResult")) {
      return(as.data.frame(obj))
    } else if (is.data.frame(obj)) {
      return(obj)
    }
    rlang::abort("Unsupported input type. Expected enrichResult, compareClusterResult, gseaResult, or data.frame.")
  }

  get_top_n <- function(df, n) {
    if (nrow(df) == 0) {
      return(df)
    }
    if ("p.adjust" %in% names(df)) {
      df <- df[order(df$p.adjust), ]
    } else if ("pvalue" %in% names(df)) {
      df <- df[order(df$pvalue), ]
    }
    utils::head(df, n)
  }

  get_genes <- function(obj, cluster = NULL) {
    if (inherits(obj, "enrichResult")) {
      return(obj@gene)
    } else if (inherits(obj, "compareClusterResult")) {
      if (!is.null(cluster) && !is.null(obj@geneClusters)) {
        if (cluster %in% names(obj@geneClusters)) {
          return(obj@geneClusters[[cluster]])
        }
      }
    }
    NULL
  }

  if (is.list(x) && !inherits(x, "enrichResult") && !inherits(x, "gseaResult") &&
    !inherits(x, "compareClusterResult") && !is.data.frame(x)) {
    if ("df" %in% names(x) && is.data.frame(x$df)) {
      return(list(Default = x))
    }
    dfs <- lapply(x, get_df)
    has_cluster <- all(vapply(dfs, function(d) "Cluster" %in% names(d), logical(1)))
    combined_df <- do.call(rbind, dfs)
    if (has_cluster) {
      df_list <- split(combined_df, combined_df$Cluster)
      res_list <- lapply(names(df_list), function(cl) {
        list(df = get_top_n(df_list[[cl]], n_pathways), genes = NULL)
      })
      names(res_list) <- names(df_list)
      return(res_list)
    } else {
      return(list(Default = list(df = get_top_n(combined_df, n_pathways), genes = NULL)))
    }
  } else {
    df <- get_df(x)
    if ("Cluster" %in% names(df)) {
      df_list <- split(df, df$Cluster)
      res_list <- lapply(names(df_list), function(cl) {
        list(df = get_top_n(df_list[[cl]], n_pathways), genes = get_genes(x, cl))
      })
      names(res_list) <- names(df_list)
      return(res_list)
    } else {
      return(list(Default = list(df = get_top_n(df, n_pathways), genes = get_genes(x))))
    }
  }
}

# ============================================================================
# PPI Helper
# ============================================================================

#' @keywords internal
.get_ppi_context_text <- function(genes, x = NULL, limit = 50) {
  if (length(genes) == 0) {
    return(NULL)
  }

  input_for_ppi <- utils::head(genes, limit)
  current_taxID <- "auto"
  if (!is.null(x) && inherits(x, "enrichResult") && !is.list(x)) {
    current_taxID <- tryCatch(getTaxID(x@organism), error = function(e) "auto")
  }

  tryCatch(
    {
      g <- getPPI(input_for_ppi, taxID = current_taxID, output = "igraph", network_type = "functional")
      if (is.null(g)) {
        return(NULL)
      }
      el <- igraph::as_data_frame(g, what = "edges")
      if (nrow(el) == 0) {
        return(NULL)
      }
      if ("score" %in% names(el)) el <- el[order(el$score, decreasing = TRUE), ]
      el_subset <- utils::head(el, limit)
      edges_text <- apply(el_subset, 1, function(row) {
        score_info <- if ("score" %in% names(row)) paste0(" (Score: ", row["score"], ")") else ""
        paste0(row["from"], " -- ", row["to"], score_info)
      })
      paste(edges_text, collapse = "\n")
    },
    error = function(e) NULL
  )
}

# ============================================================================
# Prompt Helpers (return system + user prompt pair)
# ============================================================================

.format_pathway_text <- function(df) {
  cols <- intersect(c("ID", "Description", "GeneRatio", "NES", "p.adjust", "pvalue", "geneID"), names(df))
  paste(
    apply(df[, cols, drop = FALSE], 1, function(row) {
      paste(names(row), row, sep = ": ", collapse = ", ")
    }),
    collapse = "\n"
  )
}

.build_data_sections <- function(pathway_text, context = NULL, ppi_network = NULL,
                                 fold_change = NULL, top_genes = NULL) {
  parts <- character(0)
  if (!is.null(context) && nzchar(context)) {
    parts <- c(parts, paste0("Experimental Context:\n", context))
  }
  parts <- c(parts, paste0("Top Enriched Terms:\n", pathway_text))
  if (!is.null(top_genes) && nzchar(top_genes)) {
    parts <- c(parts, paste0("Top Marker Genes (Highest Fold Change):\n", top_genes))
  }
  if (!is.null(ppi_network)) {
    parts <- c(parts, paste0("PPI Network (Edge List from STRING):\n", ppi_network))
  }
  if (!is.null(fold_change)) {
    parts <- c(parts, paste0("Top Regulated Genes (Log2 Fold Change):\n", fold_change))
  }
  paste(parts, collapse = "\n\n")
}

.interpretation_system_prompt <- function() {
  paste0(
    "You are an expert biologist and bioinformatics researcher.\n\n",
    "Analyze enrichment results using Chain-of-Thought reasoning:\n",
    "1. Source Deconvolution: Distinguish Biological Processes/Pathways (WHAT), ",
    "Upstream Regulators/TFs (WHO), and Phenotypes/Diseases (OUTCOME).\n",
    "2. Gene-Level Analysis: Identify shared key genes and unique functional modules.\n",
    "3. Causal Integration: Connect Regulators -> Processes -> Outcomes.\n",
    "4. Contextual Mapping: Map themes to the experimental context.\n",
    "5. Network Analysis: Identify functional modules and key hubs from PPI data.\n",
    "6. Network Refinement: Prune PPI to the most biologically relevant edges.\n\n",
    "GROUNDING INSTRUCTION:\n",
    "- Every claim MUST be supported by specific evidence from the provided enrichment results.\n",
    "- Cite supporting pathway names in parentheses.\n",
    "- Do not make claims that cannot be directly inferred from the provided data."
  )
}

.annotation_system_prompt <- function(cluster_id, has_prior = FALSE) {
  base <- paste0(
    "You are an expert cell biologist analyzing cluster ", cluster_id,
    " from a single-cell RNA-seq experiment.\n\n",
    "Use the following logic:\n",
    "1. Source Deconvolution: Distinguish Cell Type Markers, Biological Pathways, and Upstream TFs.\n",
    "2. Comparative Analysis (CRITICAL): Do NOT just look at the top 1 enriched term.\n",
    "   - Compare the top 3-5 enriched terms as candidates.\n",
    "   - Use specific Marker Genes to vote among candidates.\n",
    "   - Rule of Exclusion: If top term is 'Cell Type A' but gene list has markers for 'Cell Type B' ",
    "and lacks key markers for 'Cell Type A', assign 'Cell Type B'.\n",
    "3. Pathway Context: Use functional pathways to infer cell state or function.\n",
    "4. Integration: Combine marker specificity with pathway function.\n",
    "5. Network Analysis: Use PPI data to identify functional modules.\n\n",
    "Confidence Calibration:\n",
    "- High: Specific markers definitively distinguish the label.\n",
    "- Medium: Strong shared markers but weak discriminatory markers.\n",
    "- Low: Conflicting evidence.\n\n",
    "GROUNDING INSTRUCTION (STRICT):\n",
    "- NO HALLUCINATION: Do not invent markers or pathways not present in the input.\n",
    "- CITATION REQUIRED: Every conclusion must cite specific pathways or markers from the data."
  )
  if (has_prior) {
    paste0(base, "\n\nTask: Validate and refine the preliminary annotation based on the enrichment and marker evidence.")
  } else {
    paste0(base, "\n\nTask: Identify the cell type of this cluster based on the enrichment results and marker genes.")
  }
}

.phenotype_system_prompt <- function(group_id) {
  paste0(
    "You are an expert biologist characterizing the biological phenotype of group ", group_id, ".\n\n",
    "Use the following logic:\n",
    "1. Source Deconvolution: Separate observed processes (Pathways) from drivers (TFs).\n",
    "2. Synthesize enriched terms to identify the dominant biological theme.\n",
    "3. Be specific about direction or nature of the state ",
    "(e.g., 'M1 Macrophage Polarization' over 'Immune response').\n",
    "4. Assign a concise Phenotype Label.\n",
    "5. Network Analysis: Use PPI data to identify functional modules and key drivers.\n\n",
    "GROUNDING INSTRUCTION:\n",
    "- Every phenotype claim must be supported by evidence.\n",
    "- List the specific pathways or TFs that define this phenotype."
  )
}

# ============================================================================
# Core: interpret()
# ============================================================================

#' Interpret enrichment results using LLMs
#'
#' Sends enrichment results along with optional experimental context to an LLM
#' to generate a structured biological interpretation, hypothesis, and narrative
#' suitable for a publication.
#'
#' Uses `generate_object()` internally for reliable structured output with
#' automatic JSON repair, eliminating manual parsing failures.
#'
#' @param x An enrichment result object (`enrichResult`, `gseaResult`,
#'   `compareClusterResult`, or a `data.frame` with pathway columns).
#' @param context A string describing the experimental background
#'   (e.g., "scRNA-seq of mouse myocardial infarction at day 3").
#' @param n_pathways Number of top significant pathways to include. Default 20.
#' @param model Optional LLM model. When `NULL` (default), uses the aisdk
#'   package-wide default model configured via `aisdk::set_model()`.
#'   You can also supply a model ID in `provider:model` format
#'   (e.g., `"deepseek:deepseek-chat"`, `"gemini:gemini-2.5-flash"`) or a
#'   `LanguageModelV1` object. Bare model names are supported with a warning
#'   (e.g., `"deepseek-chat"`).
#' @param task Task type: "interpretation" (default), "cell_type"/"annotation",
#'   or "phenotype"/"phenotyping".
#' @param prior Optional prior knowledge or preliminary annotation to guide the task.
#' @param add_ppi Logical, whether to query STRING PPI network data. Default FALSE.
#' @param gene_fold_change Named numeric vector of log fold changes for expression context.
#' @param max_tokens Maximum tokens for the LLM response. Default 8192.
#'   Some models (especially reasoning models) may need much higher values
#'   (e.g., 16384 or more) to produce complete structured output.
#' @param temperature Sampling temperature. Default 0.3.
#' @param reduce_go Logical, whether to apply semantic-similarity GO term redundancy
#'   reduction before interpretation. Only applies to GO enrichment results.
#'   Default FALSE.
#' @param go_reduction_threshold Semantic similarity threshold (0-1).
#'   Higher values = more aggressive reduction. Default 0.7.
#' @param verbose Logical, whether to print debug messages showing raw API
#'   responses, token usage, and JSON parsing details. Default FALSE.
#'   Equivalent to setting `options(aisdk.debug = TRUE)` for the call.
#' @return An `interpretation` object (list) with task-specific fields.
#'   For "interpretation": overview, key_mechanisms, hypothesis, narrative, etc.
#'   For "annotation": cell_type, confidence, reasoning, markers, etc.
#'   For "phenotype": phenotype, confidence, reasoning, key_processes, etc.
#' @export
#' @examples
#' \dontrun{
#' # Basic usage with a data frame
#' df <- data.frame(
#'   ID = c("GO:0006915", "GO:0008284"),
#'   Description = c("apoptotic process", "positive regulation of proliferation"),
#'   GeneRatio = c("10/100", "20/100"),
#'   p.adjust = c(0.01, 0.02),
#'   geneID = c("TP53/BAX", "MYC/CCND1/CDK4")
#' )
#' res <- interpret(df,
#'   model = "deepseek:deepseek-chat",
#'   context = "Cancer proliferation study"
#' )
#' # With GO redundancy reduction
#' ego <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, ont = "BP")
#' res <- interpret(ego, reduce_go = TRUE, go_reduction_threshold = 0.7)
#' # Reuse aisdk's global default model
#' # aisdk::set_model("openai:gpt-4o-mini")
#' # res <- interpret(df, context = "Cancer proliferation study")
#' print(res)
#' }
interpret <- function(x,
                      context = NULL,
                      n_pathways = 20,
                      model = NULL,
                      task = "interpretation",
                      prior = NULL,
                      add_ppi = FALSE,
                      gene_fold_change = NULL,
                      max_tokens = 8192,
                      temperature = 0.3,
                      reduce_go = FALSE,
                      go_reduction_threshold = 0.7,
                      verbose = FALSE) {
  if (missing(x)) rlang::abort("Enrichment result 'x' is required.")

  old_debug <- getOption("aisdk.debug", FALSE)
  if (isTRUE(verbose)) {
    options(aisdk.debug = TRUE)
    on.exit(options(aisdk.debug = old_debug), add = TRUE)
  }

  # Apply rrvgo reduction if requested and input is enrichResult
  if (isTRUE(reduce_go) && is(x, "enrichResult")) {
    x <- reduce_go_redundancy(x, threshold = go_reduction_threshold)
    if (isTRUE(attr(x, "rrvgo_reduced"))) {
      message(sprintf(
        "Applied rrvgo reduction: %d -> %d GO terms (%.1f%% reduction)",
        attr(x, "rrvgo_original_n"),
        attr(x, "rrvgo_reduced_n"),
        100 * (attr(x, "rrvgo_original_n") - attr(x, "rrvgo_reduced_n")) / attr(x, "rrvgo_original_n")
      ))
    }
  }

  model <- .normalize_interpret_model(model)
  res_list <- process_enrichment_input(x, n_pathways)

  .prov <- .prov_new()

  if (length(res_list) == 0) {
    return(.attach_provenance(
      structure(
        list(overview = "No significant pathways found to interpret.", confidence = "None"),
        class = c("interpretation", "list")
      ),
      .prov_finalize(.prov, "interpret", NA_character_)
    ))
  }

  .ev_sha <- .evidence_digest_all(x, res_list)

  results <- lapply(names(res_list), function(name) {
    message(sprintf("Interpreting cluster: %s", name))
    item <- res_list[[name]]
    df <- item$df
    genes <- item$genes

    top_genes_text <- .get_top_genes_text(genes, gene_fold_change)

    if (nrow(df) == 0 && is.null(top_genes_text)) {
      res <- list(
        cluster = name,
        overview = "No significant pathways enriched and no marker genes available for interpretation.",
        confidence = "None"
      )
      class(res) <- c("interpretation", "list")
      return(res)
    }

    pathway_text <- if (nrow(df) > 0) .format_pathway_text(df) else "No significant enriched pathways found."

    current_prior <- .resolve_prior(prior, name)
    ppi_text <- .get_ppi_if_requested(add_ppi, df, genes, x)
    fc_text <- .get_fc_text(gene_fold_change, df, genes)

    user_prompt <- .build_data_sections(
      pathway_text, context, ppi_text, fc_text, top_genes_text
    )

    if (!is.null(current_prior) && nzchar(current_prior)) {
      user_prompt <- paste0(user_prompt, "\n\nPreliminary Annotation:\n", current_prior)
    }

    res <- .call_generate_object(
      model = model,
      task = task,
      cluster_id = name,
      user_prompt = user_prompt,
      has_prior = !is.null(current_prior),
      max_tokens = max_tokens,
      temperature = temperature,
      collector = .prov,
      stage = paste0("interpret:", name)
    )

    res$cluster <- name
    .postprocess_network(res)
  })

  names(results) <- names(res_list)
  .p <- .prov_finalize(.prov, "interpret", .ev_sha)

  if (length(results) == 1 && names(results)[1] == "Default") {
    return(.attach_provenance(results[[1]], .p))
  }
  class(results) <- c("interpretation_list", "list")
  .attach_provenance(results, .p)
}

# ============================================================================
# Core: interpret_agent()
# ============================================================================

## Agent 4 (Reviewer): schema + system prompt, used by interpret_agent() when review = TRUE.
## Internal helpers (not exported); audits the synthesized draft for grounding,
## overclaiming, and experimental-label fidelity.
.reviewer_schema <- function() {
  z_object(
    verdict = z_string("'pass' if every claim is grounded; 'revise' if any unsupported, overclaimed, or mislabeled claim exists"),
    grounding_score = z_string("High / Medium / Low: how well the draft's claims are supported by the provided evidence"),
    issues = z_array(
      z_object(
        claim = z_string("the specific problematic statement, quoted from the draft"),
        problem = z_string("why it is unsupported, overclaimed, or mislabeled"),
        severity = z_string("high / medium / low")
      ),
      description = "Grounding problems found; empty array if none"
    ),
    evidence_sufficiency = z_enum(
      c("sufficient", "partial", "insufficient"),
      description = "Whether the supplied evidence classes can support a claim at the resolution the draft asserts. 'insufficient' when the draft names a specific entity that the evidence cannot single out."
    ),
    supported_resolution = z_string(
      "The finest level of claim the evidence actually supports, stated plainly (e.g. 'pathway/module level: oxidative phosphorylation is perturbed' rather than 'gene level: ATP5F1C was knocked down'). Use the draft's own domain vocabulary."
    ),
    abstain_recommended = z_enum(
      c("no", "yes"),
      description = "'yes' if the specific attribution should be withheld entirely and only the supported-resolution statement reported. Independent of grounding_score: a draft can cite only real genes and still over-resolve."
    ),
    corrected_overview = z_string("Overview rewritten to remove unsupported claims, fix any mislabeled experimental groups or contrast direction, and downgrade any claim that exceeds supported_resolution to that resolution; if nothing needs fixing, faithfully restate the original"),
    adjusted_confidence = z_string("High / Medium / Low after audit"),
    warnings = z_array(z_string(), description = "Caveats a reader must be told")
  )
}

.reviewer_system_prompt <- function() {
  paste0(
    "You are 'Agent Reviewer', a skeptical bioinformatics auditor. You did NOT run the analysis; ",
    "your only job is to verify the DRAFT interpretation against the PROVIDED EVIDENCE.\n",
    "Rules:\n",
    "1. GROUNDING: every mechanism, driver, cell-type or phenotype claim must trace to an enriched pathway or gene present in the evidence. Flag any pathway or gene cited that is NOT in the evidence. Distinguish an explicitly-labelled *inferred* regulator (allowed if marked as an inference) from a gene falsely presented as observed in the data (a hallucination).\n",
    "2. CONTRAST / LABEL FIDELITY: the draft must describe the experimental groups exactly as the context states. If the context's contrast is 'A_vs_B', the draft must NOT silently relabel it (e.g. calling a susceptible-vs-control contrast 'resilient vs non-resilient').\n",
    "3. OVERCLAIMING: flag correlation-as-causation, claims of 'discovered/proven' mechanism, or phenotype attribution unsupported by enrichment.\n",
    "4. EVIDENCE SUFFICIENCY (resolution): judge whether the evidence CLASSES supplied can support a claim as specific as the draft makes, and set evidence_sufficiency, supported_resolution and abstain_recommended accordingly. This is separate from grounding: a draft may cite only genes that are present and still resolve far beyond what the evidence can single out. Guidance: enrichment of a shared process (e.g. many members of one complex or pathway) supports a claim at PATHWAY/COMPLEX level, but does not single out WHICH member was perturbed; without directionality, upstream-regulator or TF enrichment, or interaction evidence, a named single-gene cause is not supported. Naming a gene while separately declaring low confidence does NOT resolve this -- if the evidence cannot support the attribution, set abstain_recommended = 'yes' and state the supported claim instead. Do not manufacture doubt either: when the evidence does single out the entity, say sufficient.\n",
    "5. Do NOT introduce new biology, genes, or pathways. Only audit and, where needed, soften, downgrade in resolution, or correct.\n",
    "Return verdict, grounding_score, evidence_sufficiency, supported_resolution, abstain_recommended, issues, a corrected_overview, adjusted_confidence, and reader warnings."
  )
}

#' Interpret enrichment results using a multi-agent pipeline (Deep Mode)
#'
#' Employs specialized AI agents in sequence for rigorous interpretation:
#' \enumerate{
#'   \item Agent Cleaner: Filters noise and selects relevant pathways.
#'   \item Agent Detective: Identifies key regulators and functional modules.
#'   \item Agent Synthesizer: Produces a coherent biological narrative.
#'   \item Agent Reviewer (optional, \code{review = TRUE}): Audits the draft for
#'     grounding, overclaiming, and experimental-label fidelity; flags issues and
#'     returns a corrected overview.
#' }
#'
#' Uses aisdk's Agent and Session system for shared context across agents.
#'
#' @param x An enrichment result object.
#' @param context A string describing the experimental background.
#' @param n_pathways Number of top pathways to consider initially. Default 50.
#' @param model Optional LLM model. When `NULL` (default), uses the aisdk
#'   package-wide default model configured via `aisdk::set_model()`. You can
#'   also supply a model ID in `provider:model` format or a `LanguageModelV1`
#'   object. Bare model names are supported with a warning.
#' @param add_ppi Logical, whether to query PPI data. Default FALSE.
#' @param gene_fold_change Named numeric vector of log fold changes.
#' @param max_tokens Maximum tokens per agent call. Default 8192.
#' @param temperature Sampling temperature. Default 0.3.
#' @param review Logical, whether to run a final Reviewer agent (Agent 4) that
#'   audits the synthesized draft against the input evidence for grounding,
#'   overclaiming, experimental-label fidelity, and **evidence sufficiency**
#'   (whether the supplied evidence classes can support a claim as specific as
#'   the draft makes). When `TRUE` (default), the result gains a `review` field
#'   (verdict, grounding_score, evidence_sufficiency, supported_resolution,
#'   abstain_recommended, issues, adjusted_confidence, warnings) and
#'   `overview_original`; if the verdict is `"revise"`, `overview` is replaced
#'   by the corrected version. Set `FALSE` for the original three-agent
#'   behavior. Default TRUE.
#' @param verbose Logical, whether to print debug messages. Default FALSE.
#' @return An `interpretation` object with deep analysis fields plus
#'   regulatory_drivers, refined_network, and network_evidence from the
#'   detective agent. When `review = TRUE`, also includes `review` and
#'   `overview_original`. The `review` field reports grounding separately from
#'   resolution: `abstain_recommended = "yes"` flags a draft whose cited
#'   evidence is real but which resolves to a more specific entity than that
#'   evidence can single out, with `supported_resolution` giving the claim the
#'   evidence does support.
#' @export
#' @examples
#' \dontrun{
#' res <- interpret_agent(df,
#'   model = "openai:gpt-4o",
#'   context = "scRNA-seq of mouse MI day 3"
#' )
#' print(res)
#' }
interpret_agent <- function(x,
                            context = NULL,
                            n_pathways = 50,
                            model = NULL,
                            add_ppi = FALSE,
                            gene_fold_change = NULL,
                            max_tokens = 8192,
                            temperature = 0.3,
                            review = TRUE,
                            verbose = FALSE) {
  if (missing(x)) rlang::abort("Enrichment result 'x' is required.")

  old_debug <- getOption("aisdk.debug", FALSE)
  if (isTRUE(verbose)) {
    options(aisdk.debug = TRUE)
    on.exit(options(aisdk.debug = old_debug), add = TRUE)
  }

  model <- .normalize_interpret_model(model)
  res_list <- process_enrichment_input(x, n_pathways)

  .prov <- .prov_new()

  if (length(res_list) == 0) {
    return(.attach_provenance(
      structure(
        list(overview = "No significant pathways found to interpret."),
        class = c("interpretation", "list")
      ),
      .prov_finalize(.prov, "interpret_agent", NA_character_)
    ))
  }

  .ev_sha <- .evidence_digest_all(x, res_list)

  results <- lapply(names(res_list), function(name) {
    item <- res_list[[name]]
    df <- item$df
    original_genes <- item$genes
    fallback_mode <- FALSE

    if (nrow(df) == 0) {
      if (!is.null(original_genes) && length(original_genes) > 0) {
        fallback_mode <- TRUE
        rlang::warn(sprintf(
          "Cluster '%s': No enriched pathways. Falling back to gene-based interpretation.", name
        ))
        pathway_text <- paste(
          "No significant pathways enriched.",
          "Top Genes:", paste(utils::head(original_genes, 50), collapse = ", ")
        )
      } else {
        return(NULL)
      }
    } else {
      pathway_text <- .format_pathway_text(df)
    }

    ppi_text <- .get_ppi_if_requested(add_ppi, df, original_genes, x, fallback_mode)
    fc_text <- .get_fc_text(gene_fold_change, df, original_genes, fallback_mode)

    session <- ChatSession$new(model = model)

    # --- Agent 1: Cleaner ---
    cleaned_pathways <- pathway_text
    if (!fallback_mode) {
      message(sprintf("Processing cluster '%s' with Agent 1: The Cleaner...", name))
      cleaner <- Agent$new(
        name = "cleaner",
        description = "Filters noise and selects relevant pathways from enrichment results",
        system_prompt = paste0(
          "You are 'Agent Cleaner', an expert bioinformatics curator.\n",
          "Your task is to filter enriched pathways to retain only those relevant to the experimental context.\n\n",
          "Instructions:\n",
          "1. REMOVE 'housekeeping' pathways (Ribosome, Spliceosome, RNA transport) unless specifically relevant.\n",
          "2. REMOVE redundant or overly broad terms.\n",
          "3. KEEP disease-specific, tissue-specific, or phenotype-driving pathways."
        ),
        model = model
      )

      cleaner_prompt <- paste0(
        if (!is.null(context)) paste0("Context: ", context, "\n\n") else "",
        "Raw Enriched Pathways:\n", pathway_text
      )

      cleaner_res <- tryCatch(
        {
          gen <- .gen_object_stamped(
            collector = .prov, stage = "cleaner",
            model = model, prompt = cleaner_prompt,
            schema = .cleaner_schema(), schema_name = "cleaner_result",
            system = cleaner$system_prompt,
            temperature = temperature, max_tokens = max_tokens
          )
          # A 200 that does not parse is a failure, not a result. generate_object() raises no
          # condition in that case, so a bare `gen$object` hands NULL to the caller; every
          # downstream is.list()/!is.null() guard then skips silently, and for the synthesizer
          # the NULL reaches .postprocess_network(), where `class(res) <-` on NULL aborts the
          # whole run with "attempt to set an attribute on NULL" (observed 2026-08-27, both
          # review = TRUE arms of the A1 public-case ablation). Raising here routes the failure
          # to the handler below, which warns and -- for the synthesizer -- returns a marked
          # placeholder, so the call is reported rather than fatal.
          .stop_if_unparsed(gen, "Cleaner")
          gen$object
        },
        error = function(e) {
          rlang::warn(paste0("Agent Cleaner failed: ", e$message, ". Using unfiltered pathways."))
          NULL
        }
      )

      if (!is.null(cleaner_res) && !is.null(cleaner_res$kept_pathways)) {
        cleaned_pathways <- paste(
          "Selected Relevant Pathways (filtered by Agent Cleaner):",
          paste(cleaner_res$kept_pathways, collapse = ", "),
          "\nReasoning:", cleaner_res$reasoning,
          sep = "\n"
        )
        session$set_memory("cleaner_result", cleaner_res)
      }
    }

    # --- Agent 2: Detective ---
    message(sprintf("Processing cluster '%s' with Agent 2: The Detective...", name))
    detective <- Agent$new(
      name = "detective",
      description = "Identifies key regulators and functional modules from enrichment and network data",
      system_prompt = paste0(
        "You are 'Agent Detective', an expert systems biologist.\n",
        "Identify Key Drivers (Regulators) and Functional Modules from the filtered pathways and network data.\n\n",
        "Instructions:\n",
        "1. Identify potential Master Regulators (TFs, Kinases) that explain the pathways.\n",
        "2. Define Functional Modules using PPI network connections.\n",
        "3. Refine the PPI network to a core regulatory sub-network."
      ),
      model = model
    )

    detective_prompt <- paste0(
      if (!is.null(context)) paste0("Context: ", context, "\n\n") else "",
      if (fallback_mode) "WARNING: No significant enriched pathways found. Analyzing RAW GENE LISTS.\n\n" else "",
      "Filtered Pathways:\n", cleaned_pathways,
      if (!is.null(ppi_text)) paste0("\n\nPPI Network Evidence:\n", ppi_text) else "",
      if (!is.null(fc_text)) paste0("\n\nGene Fold Changes:\n", fc_text) else ""
    )

    detective_res <- tryCatch(
      {
        gen <- .gen_object_stamped(
          collector = .prov, stage = "detective",
          model = model, prompt = detective_prompt,
          schema = .detective_schema(), schema_name = "detective_result",
          system = detective$system_prompt,
          temperature = temperature, max_tokens = max_tokens
        )
        .stop_if_unparsed(gen, "Detective")
        gen$object
      },
      error = function(e) {
        rlang::warn(paste0("Agent Detective failed: ", e$message))
        NULL
      }
    )

    if (!is.null(detective_res)) {
      session$set_memory("detective_result", detective_res)
    }

    # --- Agent 3: Synthesizer ---
    message(sprintf("Processing cluster '%s' with Agent 3: The Storyteller...", name))
    detective_text <- .format_detective_report(detective_res)

    synthesizer_prompt <- paste0(
      if (!is.null(context)) paste0("Context: ", context, "\n\n") else "",
      if (fallback_mode) "WARNING: No significant enriched pathways found. Interpret with caution.\n\n" else "",
      "Data Sources:\n",
      "1. Relevant Pathways:\n", cleaned_pathways, "\n\n",
      "2. Detective's Report (Drivers & Modules):\n", detective_text
    )

    synth_system <- paste0(
      "You are 'Agent Storyteller', a senior scientific writer.\n",
      "Synthesize the findings from previous agents into a coherent biological narrative.\n\n",
      "Instructions:\n",
      "1. Write a comprehensive Overview.\n",
      "2. Explain Key Mechanisms, linking Regulators -> Modules -> Pathways.\n",
      "3. Formulate a Hypothesis.\n",
      "4. Draft a Narrative paragraph for a paper."
    )

    final_res <- tryCatch(
      {
        gen <- .gen_object_stamped(
          collector = .prov, stage = "synthesizer",
          model = model, prompt = synthesizer_prompt,
          schema = .interpretation_schema(), schema_name = "synthesis_result",
          system = synth_system,
          temperature = temperature, max_tokens = max_tokens
        )
        .stop_if_unparsed(gen, "Synthesizer")
        gen$object
      },
      error = function(e) {
        rlang::warn(paste0("Agent Synthesizer failed: ", e$message))
        list(overview = "Agent Synthesizer failed to produce structured output.", confidence = "None")
      }
    )

    if (is.list(final_res)) {
      final_res$cluster <- name
      if (fallback_mode) final_res$data_source <- "gene_list_only"
      if (!is.null(detective_res)) {
        final_res$regulatory_drivers <- detective_res$key_drivers
        final_res$refined_network <- detective_res$refined_network
        final_res$network_evidence <- detective_res$network_evidence
      }
    }

    ## Agent 4: The Reviewer -- audit the draft for grounding, overclaiming, label fidelity,
    ## and whether the evidence class present can support the resolution the draft claims.
    if (isTRUE(review) && is.list(final_res) && !is.null(final_res$overview)) {
      message(sprintf("Processing cluster '%s' with Agent 4: The Reviewer...", name))
      ## Evidence inventory: state which evidence classes are present and which are ABSENT.
      ## Without this the Reviewer cannot judge sufficiency -- it sees what was given but has no
      ## way to know what was withheld, so it audits only what is in front of it and lets a
      ## single-gene attribution built on term-level evidence pass.
      evidence_inventory <- paste0(
        "EVIDENCE INVENTORY (what the analysts did and did NOT have)\n",
        "- Enriched term table with member genes: PRESENT\n",
        "- Experimental context: ", if (!is.null(context) && nzchar(context)) "PRESENT" else "ABSENT", "\n",
        "- Gene-level directionality (log fold change): ",
        if (!is.null(fc_text)) "PRESENT" else "ABSENT (no up/down direction available)", "\n",
        "- Protein-protein interaction network: ",
        if (!is.null(ppi_text)) "PRESENT" else "ABSENT (no physical-interaction evidence)", "\n",
        "- Inferred-driver report from the Detective agent: ",
        if (!is.null(detective_text) && nzchar(detective_text))
          "PRESENT (these are INFERENCES, not observations)" else "ABSENT", "\n",
        "- Upstream regulator / TF enrichment, ChIP or perturbation evidence: ",
        "ABSENT unless it appears in the term table above\n\n"
      )
      reviewer_prompt <- paste0(
        "EVIDENCE GIVEN TO THE ANALYSTS\n",
        if (!is.null(context)) paste0("Experimental context: ", context, "\n\n") else "",
        "Enriched pathways (name and member genes):\n", cleaned_pathways, "\n\n",
        if (!is.null(detective_text) && nzchar(detective_text))
          paste0("Detective's report (inferred drivers/modules):\n", detective_text, "\n\n") else "",
        evidence_inventory,
        "DRAFT INTERPRETATION TO AUDIT\n",
        "Overview: ", final_res$overview %||% "", "\n",
        "Key mechanisms: ", final_res$key_mechanisms %||% "", "\n",
        "Hypothesis: ", final_res$hypothesis %||% "", "\n",
        "Narrative: ", final_res$narrative %||% "", "\n",
        "Regulatory drivers: ", paste(unlist(final_res$regulatory_drivers), collapse = ", "), "\n\n",
        "Audit the draft strictly against the evidence, and judge separately whether the ",
        "evidence classes listed in the inventory can support the resolution the draft claims. ",
        "A draft that states LOW confidence and then asserts a specific attribution anyway is ",
        "exactly the case rule 4 exists to catch: report it, set the resolution the evidence ",
        "does support, and rewrite the overview at that resolution."
      )
      review_res <- tryCatch(
        {
          gen <- .gen_object_stamped(
            collector = .prov, stage = "reviewer",
            model = model, prompt = reviewer_prompt,
            schema = .reviewer_schema(), schema_name = "reviewer_result",
            system = .reviewer_system_prompt(),
            temperature = temperature, max_tokens = max_tokens
          )
          .stop_if_unparsed(gen, "Reviewer")
          gen$object
        },
        error = function(e) {
          rlang::warn(paste0("Agent Reviewer failed: ", e$message))
          NULL
        }
      )
      if (!is.null(review_res)) {
        final_res$overview_original <- final_res$overview
        final_res$review <- review_res
        if (!is.null(review_res$verdict) && tolower(review_res$verdict) == "revise" &&
            !is.null(review_res$corrected_overview) && nzchar(review_res$corrected_overview)) {
          final_res$overview <- review_res$corrected_overview
        }
      }
    }

    .postprocess_network(final_res)
  })

  results <- Filter(Negate(is.null), results)
  names(results) <- vapply(results, function(r) r$cluster %||% "Unknown", character(1))
  .p <- .prov_finalize(.prov, "interpret_agent", .ev_sha,
                       extra = list(review_enabled = isTRUE(review)))

  if (length(results) == 1 && names(results)[1] == "Default") {
    return(.attach_provenance(results[[1]], .p))
  }
  class(results) <- c("interpretation_list", "list")
  .attach_provenance(results, .p)
}

# ============================================================================
# Core: interpret_hierarchical()
# ============================================================================

#' Interpret enrichment results using a hierarchical strategy
#'
#' First interprets major clusters to establish lineage context, then interprets
#' sub-clusters with hierarchical constraints from the major cluster annotations.
#'
#' @param x_minor Enrichment result for sub-clusters.
#' @param x_major Enrichment result for major clusters.
#' @param mapping A named vector mapping sub-cluster IDs to major cluster IDs.
#' @param model Optional LLM model. When `NULL` (default), uses the aisdk
#'   package-wide default model configured via `aisdk::set_model()`. You can
#'   also supply a model ID in `provider:model` format or a `LanguageModelV1`
#'   object. Bare model names are supported with a warning.
#' @param task Task type, default "cell_type".
#' @param max_tokens Maximum tokens. Default 8192.
#' @param temperature Sampling temperature. Default 0.3.
#' @return An `interpretation_list` object.
#' @export
interpret_hierarchical <- function(x_minor,
                                   x_major,
                                   mapping,
                                   model = NULL,
                                   task = "cell_type",
                                   max_tokens = 8192,
                                   temperature = 0.3) {
  model <- .normalize_interpret_model(model)

  message("Step 1: Interpreting Major Clusters to establish lineage context...")
  res_major <- interpret(
    x_major,
    context = NULL, model = model, task = "cell_type",
    max_tokens = max_tokens, temperature = temperature
  )

  message("Step 2: Interpreting Sub-clusters using hierarchical constraints...")
  res_list_minor <- process_enrichment_input(x_minor, n_pathways = 20)

  results <- lapply(names(res_list_minor), function(name) {
    specific_context <- NULL
    if (name %in% names(mapping)) {
      major_id <- mapping[[name]]
      major_info <- .extract_major_info(res_major, major_id)
      if (!is.null(major_info) && !is.null(major_info$cell_type)) {
        specific_context <- paste0(
          "Hierarchical Constraint: This cluster is a confirmed subcluster of the '",
          major_info$cell_type,
          "' lineage. Focus on distinguishing the specific subtype or state within this lineage."
        )
      }
    }

    if (is.null(specific_context)) {
      rlang::warn(paste("No major lineage context found for sub-cluster:", name))
    }

    res <- interpret(
      res_list_minor[[name]],
      context = specific_context, model = model,
      task = task, max_tokens = max_tokens, temperature = temperature
    )
    if (is.list(res)) res$cluster <- name
    res
  })

  names(results) <- names(res_list_minor)
  class(results) <- c("interpretation_list", "list")
  results
}

# ============================================================================
# Internal Helpers
# ============================================================================

.diagnose_failure <- function(finish_reason, raw_text, max_tokens) {
  raw_len <- nchar(raw_text %||% "")

  if (finish_reason == "length" && raw_len == 0) {
    return(list(
      message = paste0(
        "Token limit exhausted (max_tokens=", max_tokens, "). ",
        "The model used all tokens on internal reasoning and produced no output."
      ),
      suggestion = paste0(
        "Increase max_tokens. Try: max_tokens = ", max(max_tokens * 4, 16384)
      ),
      fallback_overview = paste0(
        "Model ran out of tokens (max_tokens=", max_tokens,
        "). Re-run with a higher max_tokens value."
      )
    ))
  }

  if (finish_reason == "length" && raw_len > 0) {
    return(list(
      message = paste0(
        "Token limit exhausted (max_tokens=", max_tokens, "). ",
        "The JSON output was truncated, causing parse failure."
      ),
      suggestion = paste0(
        "Increase max_tokens. Try: max_tokens = ", max(max_tokens * 2, 8192)
      ),
      fallback_overview = raw_text
    ))
  }

  if (raw_len > 0) {
    return(list(
      message = "Model returned text but it could not be parsed as valid JSON.",
      suggestion = "The model may not follow JSON instructions well. Try a different model or lower temperature.",
      fallback_overview = raw_text
    ))
  }

  list(
    message = "Model returned an empty response.",
    suggestion = "Try a different model, or check that your API key has access to this model.",
    fallback_overview = "Model returned empty response. Try a different model or increase max_tokens."
  )
}

.get_top_genes_text <- function(genes, gene_fold_change) {
  if (is.null(genes) || length(genes) == 0) {
    return(NULL)
  }
  if (!is.null(gene_fold_change)) {
    common <- intersect(genes, names(gene_fold_change))
    if (length(common) > 0) {
      fc <- gene_fold_change[common]
      top_up <- utils::head(names(fc[order(fc, decreasing = TRUE)]), 20)
      return(paste(top_up, collapse = ", "))
    }
  }
  paste(utils::head(genes, 20), collapse = ", ")
}

.resolve_prior <- function(prior, name) {
  if (is.null(prior)) {
    return(NULL)
  }
  if (length(prior) == 1 && is.null(names(prior))) {
    return(prior)
  }
  if (name %in% names(prior)) {
    return(prior[[name]])
  }
  NULL
}

.get_ppi_if_requested <- function(add_ppi, df, genes, x, fallback_mode = FALSE) {
  if (!add_ppi) {
    return(NULL)
  }
  if (fallback_mode) {
    all_genes <- genes
  } else {
    all_genes <- unique(unlist(strsplit(as.character(df$geneID), "/")))
    if (length(all_genes) == 0 && !is.null(genes)) all_genes <- utils::head(genes, 50)
  }
  if (length(all_genes) > 0) .get_ppi_context_text(all_genes, x) else NULL
}

.get_fc_text <- function(gene_fold_change, df, genes, fallback_mode = FALSE) {
  if (is.null(gene_fold_change)) {
    return(NULL)
  }
  if (fallback_mode) {
    all_genes <- genes
  } else {
    all_genes <- unique(unlist(strsplit(as.character(df$geneID), "/")))
    if (length(all_genes) == 0 && !is.null(genes)) all_genes <- genes
  }
  common_genes <- intersect(all_genes, names(gene_fold_change))
  if (length(common_genes) == 0) {
    return(NULL)
  }
  fc_subset <- gene_fold_change[common_genes]
  fc_subset <- fc_subset[order(abs(fc_subset), decreasing = TRUE)]
  top_fc <- utils::head(fc_subset, 20)
  paste(names(top_fc), round(top_fc, 2), sep = ":", collapse = ", ")
}

# ============================================================================
# Execution provenance
# ============================================================================
#
# Why this exists. A remote model is not a versionable artifact the way a Bioconductor
# dependency is: the string "gpt-5.5" names whatever the gateway serves today. Two failures
# in this project made the consequence concrete. (1) A decoding parameter the caller set was
# silently dropped before reaching the model, and nothing in the returned object would have
# revealed it. (2) A batch of runs recorded an HTTP 401 error string as though it were the
# interpretation, because the call raised no condition the caller could catch and the result
# object carried no status. Both are invisible to an analyst holding only the result.
#
# The stamp below closes that gap: every model call in this module records what was asked
# for, what came back, and what evidence it was asked about, hashed so that a claim can be
# tied to the exact input that licensed it.

#' @keywords internal
.sha256 <- function(x) {
  if (is.null(x)) return(NA_character_)
  s <- paste(as.character(x), collapse = "\n")
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(s, algo = "sha256", serialize = FALSE))
  }
  NA_character_
}

#' @keywords internal
.model_identity <- function(model) {
  out <- list(model_id = NA_character_, provider = NA_character_, endpoint = NA_character_)
  if (is.null(model)) return(out)
  if (is.character(model) && length(model) == 1L) {
    out$model_id <- model
    parts <- strsplit(model, ":", fixed = TRUE)[[1]]
    if (length(parts) > 1L) {
      out$provider <- parts[1]
      out$model_id <- paste(parts[-1], collapse = ":")
    }
    return(out)
  }
  # aisdk LanguageModelV1 (or a custom provider's model object): introspect defensively,
  # since the field names are not part of a stable public contract.
  grab <- function(...) {
    for (f in list(...)) {
      v <- tryCatch(f(), error = function(e) NULL)
      if (!is.null(v) && is.character(v) && length(v) == 1L && nzchar(v)) return(v)
    }
    NA_character_
  }
  out$model_id <- grab(function() model$modelId, function() model$model_id,
                       function() model$model, function() attr(model, "modelId"))
  out$provider <- grab(function() model$provider, function() model$providerName,
                       function() model$provider_name, function() attr(model, "provider"))
  url <- grab(function() model$config$baseURL, function() model$config$base_url,
              function() model$baseURL, function() model$base_url)
  if (!is.na(url)) {
    # host only: the full URL can carry a key in a query string
    out$endpoint <- sub("^([a-z]+://[^/]+).*$", "\\1", url)
  }
  out
}

#' @keywords internal
.evidence_digest <- function(x, df = NULL, genes = NULL) {
  # The evidence set is what makes grounding a set-membership test rather than an opinion.
  # Hashing it lets a reader prove that a given interpretation was produced from a given
  # enrichment table, and detect the case where the object was edited in between.
  bits <- character(0)
  if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
    for (col in c("ID", "Description", "p.adjust", "geneID", "core_enrichment")) {
      if (col %in% names(df)) bits <- c(bits, paste(df[[col]], collapse = "|"))
    }
    bits <- c(bits, paste0("nrow=", nrow(df)))
  }
  if (!is.null(genes) && length(genes) > 0) {
    bits <- c(bits, paste(sort(unique(as.character(genes))), collapse = "|"))
  }
  if (length(bits) == 0) return(NA_character_)
  .sha256(bits)
}

#' @keywords internal
.prov_new <- function() {
  e <- new.env(parent = emptyenv())
  e$calls <- list()
  e$t_start <- Sys.time()
  e
}

#' Record one model call. Never throws: provenance must not be able to break a run.
#' @keywords internal
.prov_record <- function(collector, stage, model, system, prompt, schema_name,
                         temperature, max_tokens, gen = NULL, t0 = NULL,
                         error_message = NULL) {
  if (is.null(collector)) return(invisible(NULL))
  tryCatch({
    id <- .model_identity(model)
    raw <- tryCatch(gen$raw_text, error = function(e) NULL)
    obj <- tryCatch(gen$object, error = function(e) NULL)
    fin <- tryCatch(gen$finish_reason, error = function(e) NULL)
    usage <- tryCatch(gen$usage, error = function(e) NULL)
    # The model the caller ASKED for and the model that ANSWERED are different facts. Most
    # of this project's runs went through a relay whose catalogue changed under a stable
    # interface, so the requested identifier is not evidence of what served the request.
    # Where the provider echoes a model id in the response body, record it; where it does
    # not, leave NA rather than copying the request over it.
    served <- tryCatch(gen$raw_response$model, error = function(e) NULL)
    if (is.null(served)) served <- tryCatch(gen$raw_response$body$model, error = function(e) NULL)
    if (!is.character(served) || length(served) != 1L) served <- NA_character_
    status <- if (!is.null(error_message)) "call_failed" else if (!is.null(obj)) "ok" else "unparsed"
    collector$calls[[length(collector$calls) + 1L]] <- list(
      stage = stage,
      status = status,
      error_message = error_message %||% NA_character_,
      model_requested = id$model_id,
      model_served = served,
      model_matches_request = if (is.na(served)) NA else identical(served, id$model_id),
      provider = id$provider,
      endpoint = id$endpoint,
      schema_name = schema_name,
      # what was asked for
      temperature_requested = temperature,
      max_tokens_requested = max_tokens,
      # what came back. The absence of an echo of the requested parameters is itself the
      # finding: aisdk/gateway responses do not report the settings actually applied, so
      # `*_applied` stays NA and the caller can see that it is unverified rather than assume
      # the request was honoured.
      temperature_applied = NA_real_,
      max_tokens_applied = NA_real_,
      finish_reason = fin %||% NA_character_,
      response_chars = if (is.null(raw)) NA_integer_ else nchar(raw),
      prompt_tokens = usage$prompt_tokens %||% NA_integer_,
      completion_tokens = usage$completion_tokens %||% NA_integer_,
      total_tokens = usage$total_tokens %||% NA_integer_,
      system_sha256 = .sha256(system),
      prompt_sha256 = .sha256(prompt),
      response_sha256 = .sha256(raw),
      runtime_sec = if (is.null(t0)) NA_real_ else
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
      timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
  }, error = function(e) invisible(NULL))
  invisible(NULL)
}

#' @keywords internal
.pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                 error = function(e) NA_character_)

#' @keywords internal
.prov_finalize <- function(collector, entry_point, evidence_sha256 = NA_character_,
                           extra = list()) {
  if (is.null(collector)) return(NULL)
  calls <- collector$calls
  st <- vapply(calls, function(c) c$status, character(1))
  prov <- c(list(
    entry_point = entry_point,
    ok = length(st) > 0 && all(st == "ok"),
    n_calls = length(calls),
    n_calls_failed = sum(st != "ok"),
    evidence_sha256 = evidence_sha256,
    started_utc = format(collector$t_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    runtime_sec = round(as.numeric(difftime(Sys.time(), collector$t_start, units = "secs")), 2),
    software = list(
      r = paste0(R.version$major, ".", R.version$minor),
      clusterProfiler = .pkg_ver("clusterProfiler"),
      aisdk = .pkg_ver("aisdk"),
      DOSE = .pkg_ver("DOSE"),
      GOSemSim = .pkg_ver("GOSemSim")
    ),
    calls = calls
  ), extra)
  class(prov) <- c("interpret_provenance", "list")
  prov
}

#' Execution provenance of an interpretation
#'
#' Returns the execution record attached to a result of [interpret()],
#' [interpret_agent()] or [interpret_hierarchical()]: the model actually addressed, the
#' decoding parameters requested (and whether the response confirmed them), SHA-256 digests
#' of the system prompt, user prompt, raw response and input evidence, token usage,
#' `finish_reason`, per-call status and package versions.
#'
#' A remote model is not a versionable dependency. Without this record an interpretation
#' cannot be reproduced, and a call that failed at the transport layer is indistinguishable
#' from a model that had nothing to say -- both arrive as an ordinary result object. Check
#' `provenance(x)$ok` before treating any interpretation as data.
#'
#' @param x An `interpretation`, `interpretation_list` or agent result.
#' @return An `interpret_provenance` list, or `NULL` if the object carries no stamp.
#' @export
provenance <- function(x) {
  p <- attr(x, "provenance")
  if (!is.null(p)) return(p)
  if (is.list(x) && !is.null(x$provenance)) return(x$provenance)
  NULL
}

#' @export
print.interpret_provenance <- function(x, ...) {
  cat("<interpret provenance>\n")
  cat(sprintf("  entry point : %s\n", x$entry_point))
  cat(sprintf("  status      : %s (%d call%s, %d failed)\n",
              if (isTRUE(x$ok)) "ok" else "INCOMPLETE",
              x$n_calls, if (x$n_calls == 1) "" else "s", x$n_calls_failed))
  cat(sprintf("  started     : %s  (%.1fs)\n", x$started_utc, x$runtime_sec))
  cat(sprintf("  evidence    : sha256 %s\n", substr(x$evidence_sha256, 1, 16)))
  for (c in x$calls) {
    served <- if (is.na(c$model_served)) "not reported" else
      if (isTRUE(c$model_matches_request)) "as requested" else paste0("SERVED ", c$model_served)
    cat(sprintf("  - %-11s %-10s req=%s%s (%s)  T=%s max_tokens=%s  finish=%s  %sc  %.1fs\n",
                c$stage, c$status,
                if (is.na(c$provider)) "" else paste0(c$provider, ":"),
                c$model_requested, served,
                format(c$temperature_requested), format(c$max_tokens_requested),
                c$finish_reason, format(c$response_chars), c$runtime_sec))
    if (!is.na(c$error_message)) cat(sprintf("      error: %s\n", c$error_message))
  }
  cat(sprintf("  software    : R %s / clusterProfiler %s / aisdk %s\n",
              x$software$r, x$software$clusterProfiler, x$software$aisdk))
  cat("  note        : *_applied are NA because the response does not echo the settings\n")
  cat("                actually used; requested != verified. A model id shown as\n")
  cat("                'not reported' means the provider returned none -- the requested\n")
  cat("                identifier is then a request, not a record of what answered.\n")
  invisible(x)
}

#' Digest of the whole evidence set handed to one entry-point call.
#' @keywords internal
.evidence_digest_all <- function(x, res_list) {
  bits <- character(0)
  for (nm in names(res_list)) {
    it <- res_list[[nm]]
    d <- .evidence_digest(x, it$df, it$genes)
    bits <- c(bits, paste0(nm, "=", if (is.na(d)) "" else d))
  }
  if (length(bits) == 0) return(NA_character_)
  .sha256(bits)
}

#' @keywords internal
.attach_provenance <- function(obj, prov) {
  if (is.null(prov)) return(obj)
  attr(obj, "provenance") <- prov
  obj
}

#' @keywords internal
.gen_object_stamped <- function(collector, stage, model, prompt, schema, schema_name,
                                system, temperature, max_tokens) {
  t0 <- Sys.time()
  # Capture through the tryCatch return value rather than by assignment from inside it.
  # `gen <<- ...` here would be a bug: `<<-` starts its search in the PARENT environment, so
  # it skips the local binding, silently creates a global, and leaves this function returning
  # NULL on every call -- which reaches the caller as a model that produced nothing.
  cap <- tryCatch(
    list(gen = generate_object(model = model, prompt = prompt, schema = schema,
                               schema_name = schema_name, system = system,
                               temperature = temperature, max_tokens = max_tokens),
         err = NULL),
    error = function(e) list(gen = NULL, err = conditionMessage(e))
  )
  gen <- cap$gen
  err <- cap$err
  .prov_record(collector, stage, model, system, prompt, schema_name,
               temperature, max_tokens, gen = gen, t0 = t0, error_message = err)
  # Re-raise so each call site keeps the error handling it already had; the stamp is
  # recorded first, so a failed call is still visible in provenance().
  if (!is.null(err)) stop(err, call. = FALSE)
  gen
}

.call_generate_object <- function(model, task, cluster_id, user_prompt,
                                  has_prior = FALSE, max_tokens = 8192,
                                  temperature = 0.3, collector = NULL,
                                  stage = "interpret") {
  if (task %in% c("annotation", "cell_type")) {
    sys <- .annotation_system_prompt(cluster_id, has_prior)
    schema <- if (has_prior) .annotation_refinement_schema() else .annotation_schema()
    schema_name <- "annotation_result"
  } else if (task %in% c("phenotype", "phenotyping")) {
    sys <- .phenotype_system_prompt(cluster_id)
    schema <- .phenotype_schema()
    schema_name <- "phenotype_result"
  } else {
    sys <- .interpretation_system_prompt()
    schema <- .interpretation_schema()
    schema_name <- "interpretation_result"
  }

  debug <- isTRUE(getOption("aisdk.debug", FALSE))

  # `model` may be a character id OR an aisdk LanguageModelV1, which is an ENVIRONMENT. Pasting an
  # environment into a message raises "cannot coerce type 'environment' to vector of type
  # 'character'". Every place this function names the model in human-readable text must therefore
  # go through .model_identity(), never through `model` directly. The two sites in the FAILURE
  # branches below are the ones that mattered: they turned every failed call into a bogus
  # input-contract error and destroyed the real diagnosis with it.
  model_label <- {
    id <- tryCatch(.model_identity(model)$model_id, error = function(e) NA_character_)
    if (is.null(id) || is.na(id)) "<unnamed model object>" else id
  }

  if (debug) {
    message(
      "[DEBUG] .call_generate_object: model=", model_label,
      " task=", task, " cluster=", cluster_id, " max_tokens=", max_tokens
    )
  }

  result <- tryCatch(
    {
      gen <- .gen_object_stamped(
        collector = collector, stage = stage,
        model = model, prompt = user_prompt, schema = schema,
        schema_name = schema_name, system = sys,
        temperature = temperature, max_tokens = max_tokens
      )

      if (debug) {
        message("[DEBUG] generate_object result:")
        message("[DEBUG]   finish_reason: ", gen$finish_reason %||% "NULL")
        message("[DEBUG]   raw_text length: ", nchar(gen$raw_text %||% ""), " chars")
        message("[DEBUG]   object is NULL: ", is.null(gen$object))
        if (!is.null(gen$usage)) {
          message(
            "[DEBUG]   usage: prompt=", gen$usage$prompt_tokens,
            " completion=", gen$usage$completion_tokens,
            " total=", gen$usage$total_tokens
          )
        }
        if (!is.null(gen$raw_text) && nchar(gen$raw_text) > 0) {
          preview <- substr(gen$raw_text, 1, min(500, nchar(gen$raw_text)))
          message(
            "[DEBUG]   raw_text preview:\n", preview,
            if (nchar(gen$raw_text) > 500) "\n... [truncated]" else ""
          )
        }
      }

      if (!is.null(gen$object)) {
        obj <- gen$object
        class(obj) <- c("interpretation", class(obj))
        obj
      } else {
        raw <- gen$raw_text %||% ""
        finish <- gen$finish_reason %||% "unknown"
        raw_preview <- if (nchar(raw) > 200) paste0(substr(raw, 1, 200), "...") else raw

        diagnosis <- .diagnose_failure(finish, raw, max_tokens)

        warn_parts <- c(
          paste0("generate_object() returned NULL for cluster '", cluster_id, "'."),
          "i" = paste0("Model: ", model_label),
          "i" = paste0("finish_reason: ", finish),
          "i" = paste0(
            "raw_text (", nchar(raw), " chars): ",
            if (nzchar(raw_preview)) raw_preview else "<empty>"
          )
        )
        warn_parts <- c(warn_parts, "!" = diagnosis$message)
        if (!is.null(diagnosis$suggestion)) {
          warn_parts <- c(warn_parts, ">" = diagnosis$suggestion)
        }
        rlang::warn(warn_parts)

        fallback_text <- if (nzchar(raw)) raw else diagnosis$fallback_overview
        list(
          overview = fallback_text,
          confidence = "Low",
          cluster = cluster_id
        )
      }
    },
    error = function(e) {
      rlang::warn(c(
        paste0("LLM call failed for cluster '", cluster_id, "': ", e$message),
        "i" = paste0("Model: ", model_label),
        "i" = "Tip: Re-run with verbose=TRUE for full debug output"
      ))
      res <- list(
        overview = paste0("LLM call failed: ", e$message),
        confidence = "None",
        cluster = cluster_id
      )
      class(res) <- c("interpretation", "list")
      res
    }
  )

  result
}

.format_detective_report <- function(detective_res) {
  if (is.null(detective_res) || !is.list(detective_res)) {
    return("No detective report available.")
  }
  key_drivers <- if (!is.null(detective_res$key_drivers)) {
    paste(detective_res$key_drivers, collapse = ", ")
  } else {
    "None identified"
  }
  modules <- if (!is.null(detective_res$functional_modules)) {
    paste(detective_res$functional_modules, collapse = ", ")
  } else {
    "None identified"
  }
  evidence <- detective_res$network_evidence %||% "None provided"
  paste0(
    "Key Drivers: ", key_drivers, "\n",
    "Functional Modules: ", modules, "\n",
    "Network Evidence: ", evidence
  )
}

#' Raise when a model call returned successfully but produced no parsable object.
#'
#' `generate_object()` signals nothing when the transport succeeds and the body fails to
#' validate against the schema: it returns normally with `object = NULL` and the raw text
#' intact. Every agent in `interpret_agent()` therefore has to decide for itself whether a
#' missing object is "this agent had nothing to add" or "this call failed", and until
#' 2026-08-27 all four treated it as the former. `finish_reason` is carried into the message
#' because it is the one field that distinguishes a truncated response from a well-formed
#' answer the schema rejected.
#' @keywords internal
.stop_if_unparsed <- function(gen, agent) {
  if (!is.null(gen) && !is.null(gen$object)) return(invisible(TRUE))
  fr <- tryCatch(gen$finish_reason, error = function(e) NULL) %||% "NA"
  n  <- tryCatch(nchar(gen$raw_text %||% ""), error = function(e) 0L)
  stop(sprintf("Agent %s returned no parsable object (finish_reason=%s, raw_text=%d chars)",
               agent, fr, n), call. = FALSE)
}

.postprocess_network <- function(res) {
  # Defence in depth for the failure above. `class(res) <- ...` on NULL raises "attempt to set
  # an attribute on NULL" and aborts the entire run, so a cluster whose synthesizer produced
  # nothing took the whole call down with it instead of being reported as an empty result.
  # Every caller of this function is at the end of a pipeline that may legitimately have
  # produced nothing; none of them should have to crash to say so.
  if (is.null(res)) {
    res <- list(overview = "No interpretation was produced for this unit.", confidence = "None")
  }
  if (!is.list(res) || is.null(res$refined_network)) {
    class(res) <- union("interpretation", class(res))
    return(res)
  }
  rn <- res$refined_network
  rn_df <- tryCatch(
    {
      if (is.data.frame(rn)) {
        rn
      } else if (is.list(rn) && length(rn) > 0 && all(vapply(rn, is.list, logical(1)))) {
        do.call(rbind, lapply(rn, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )

  if (!is.null(rn_df) && nrow(rn_df) > 0) {
    colnames(rn_df)[colnames(rn_df) == "source"] <- "from"
    colnames(rn_df)[colnames(rn_df) == "target"] <- "to"
    if ("from" %in% names(rn_df) && "to" %in% names(rn_df)) {
      if (requireNamespace("igraph", quietly = TRUE)) {
        res$network <- igraph::graph_from_data_frame(rn_df, directed = FALSE)
      }
    }
  }
  class(res) <- union("interpretation", class(res))
  res
}

.extract_major_info <- function(res_major, major_id) {
  if (inherits(res_major, "interpretation_list") && major_id %in% names(res_major)) {
    return(res_major[[major_id]])
  }
  if (inherits(res_major, "interpretation")) {
    if (!is.null(res_major$cluster) && res_major$cluster == major_id) {
      return(res_major)
    }
    if (is.null(res_major$cluster)) {
      return(res_major)
    }
  }
  NULL
}

# ============================================================================
# Print Methods
# ============================================================================

#' @method print interpretation
#' @export
print.interpretation <- function(x, ...) {
  if (!is.null(x$cell_type)) {
    cat("## Cell Type Annotation\n\n")
    if (!is.null(x$cluster)) cat(sprintf("### Cluster: %s\n\n", x$cluster))
    cat(sprintf("**Cell Type:** %s\n", x$cell_type))
    if (!is.null(x$refinement_status)) cat(sprintf("**Status:** %s\n", x$refinement_status))
    cat(sprintf("**Confidence:** %s\n", x$confidence))
    cat("\n**Reasoning:**\n", x$reasoning, "\n")
    if (!is.null(x$markers)) {
      cat("\n**Supporting Markers/Pathways:**\n")
      cat(paste("-", unlist(x$markers), collapse = "\n"), "\n")
    }
    cat("\n")
    return(invisible(x))
  }

  if (!is.null(x$phenotype)) {
    cat("## Phenotype Characterization\n\n")
    if (!is.null(x$cluster)) cat(sprintf("### Group/Cluster: %s\n\n", x$cluster))
    cat(sprintf("**Phenotype:** %s\n", x$phenotype))
    cat(sprintf("**Confidence:** %s\n", x$confidence))
    cat("\n**Reasoning:**\n", x$reasoning, "\n")
    if (!is.null(x$key_processes)) {
      cat("\n**Key Processes:**\n")
      cat(paste("-", unlist(x$key_processes), collapse = "\n"), "\n")
    }
    cat("\n")
    return(invisible(x))
  }

  cat("## Interpretation Result\n\n")
  if (!is.null(x$cluster)) cat(sprintf("### Cluster: %s\n\n", x$cluster))

  if (!is.null(x$overview)) {
    cat("### 1. Overview\n")
    cat(x$overview, "\n\n")
  }

  if (!is.null(x$regulatory_drivers)) {
    cat("### 2. Regulatory Drivers\n")
    drivers <- unlist(x$regulatory_drivers)
    if (length(drivers) > 1) {
      cat(paste("-", drivers, collapse = "\n"), "\n\n")
    } else {
      cat(drivers, "\n\n")
    }
  }

  if (!is.null(x$key_mechanisms)) {
    cat("### 3. Key Mechanisms\n")
    if (is.list(x$key_mechanisms) && !is.null(names(x$key_mechanisms))) {
      for (mname in names(x$key_mechanisms)) {
        cat(sprintf("#### %s\n", mname))
        mechanism <- x$key_mechanisms[[mname]]
        if (is.list(mechanism)) {
          if (!is.null(mechanism$explanation)) cat(mechanism$explanation, "\n")
          if (!is.null(mechanism$pathways)) {
            cat("**Pathways:** ", paste(mechanism$pathways, collapse = ", "), "\n")
          }
          if (!is.null(mechanism$genes)) {
            cat(
              "**Key Genes:** ", paste(utils::head(mechanism$genes, 10), collapse = ", "),
              if (length(mechanism$genes) > 10) "..." else "", "\n"
            )
          }
        } else {
          cat(mechanism, "\n")
        }
        cat("\n")
      }
    } else {
      cat(as.character(x$key_mechanisms), "\n\n")
    }
  }

  if (!is.null(x$crosstalk)) {
    cat("### 4. Crosstalk & Interactions\n")
    cat(x$crosstalk, "\n\n")
  }

  if (!is.null(x$hypothesis)) {
    cat("### 5. Hypothesis\n")
    if (is.list(x$hypothesis)) {
      if (!is.null(x$hypothesis$what)) cat("**Observation (What):** ", x$hypothesis$what, "\n\n")
      if (!is.null(x$hypothesis$so_what)) cat("**Implication (So What):** ", x$hypothesis$so_what, "\n\n")
    } else {
      cat(x$hypothesis, "\n\n")
    }
  }

  if (!is.null(x$narrative)) {
    cat("### 6. Narrative Draft\n")
    cat(x$narrative, "\n\n")
  }

  if (!is.null(x$network) && requireNamespace("igraph", quietly = TRUE) && inherits(x$network, "igraph")) {
    cat("### 7. Refined Regulatory Network\n")
    el <- igraph::as_data_frame(x$network, what = "edges")
    if (nrow(el) > 0) {
      cat("Key Interactions:\n")
      for (i in seq_len(nrow(el))) {
        interaction_type <- if ("interaction" %in% names(el)) paste0(" (", el[i, "interaction"], ")") else ""
        reason <- if ("reason" %in% names(el)) paste0(" - ", el[i, "reason"]) else ""
        cat(sprintf("  %s -- %s%s%s\n", el[i, "from"], el[i, "to"], interaction_type, reason))
      }
      cat("\n")
    }
    if (!is.null(x$network_evidence)) {
      cat("**Network Evidence:**\n")
      cat(x$network_evidence, "\n\n")
    }
  }

  invisible(x)
}

#' @method print interpretation_list
#' @export
print.interpretation_list <- function(x, ...) {
  cat("# Enrichment Interpretation / Annotation Report\n\n")
  for (i in seq_along(x)) {
    print(x[[i]])
    cat("---\n\n")
  }
  invisible(x)
}

# ============================================================================
# rrvgo Integration: Optional GO Term Redundancy Reduction
# ============================================================================

#' Reduce GO term redundancy by semantic similarity
#'
#' Applies semantic similarity-based clustering to remove redundant GO terms
#' before LLM interpretation. This algorithmic preprocessing complements the
#' multi-agent interpretation layer by reducing input noise and token count.
#'
#' Implemented directly on GOSemSim, which is already an Imports dependency.
#' The `rrvgo_*` attribute names are kept for backwards compatibility with
#' existing ablation outputs and predate the switch away from rrvgo.
#'
#' @param enrich_result An enrichResult object from enrichGO()
#' @param threshold Semantic similarity threshold (0-1). Higher values are more
#'   aggressive. Default 0.7 means terms with >70% similarity are clustered.
#' @param orgdb OrgDb package name (e.g., "org.Hs.eg.db"). If NULL, auto-detected
#'   from the enrichResult organism field.
#' @param method Similarity measure passed to [GOSemSim::goSim()]. Default "Wang",
#'   which reads the GO DAG directly and needs no corpus statistics. The
#'   information-content measures "Resnik", "Lin", "Rel" and "Jiang" are accepted
#'   but require `godata(computeIC = TRUE)` and are substantially slower.
#' @return A reduced enrichResult object with redundant terms removed. The object
#'   carries attributes recording the reduction: `rrvgo_reduced`, `rrvgo_original_n`,
#'   `rrvgo_reduced_n`, `rrvgo_threshold`, `rrvgo_method`.
#' @keywords internal
#' @examples
#' \dontrun{
#' ego <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, ont = "BP")
#' ego_reduced <- reduce_go_redundancy(ego, threshold = 0.7)
#' message(sprintf("Reduced from %d to %d terms",
#'   attr(ego_reduced, "rrvgo_original_n"),
#'   attr(ego_reduced, "rrvgo_reduced_n")))
#' }
reduce_go_redundancy <- function(enrich_result,
                                 threshold = 0.7,
                                 orgdb = NULL,
                                 method = "Wang") {
  # Semantic-similarity redundancy reduction over an enrichGO result.
  #
  # Implemented directly on GOSemSim rather than on rrvgo. rrvgo is the obvious dependency
  # and we tried it first, but installing it on the execution node failed: its chain runs
  # rrvgo -> umap -> reticulate, and the node cannot complete it. GOSemSim -- which rrvgo
  # itself sits on -- is already present, and supplies everything this step needs.
  #
  # `method` defaults to Wang deliberately. The information-content measures (Rel, Resnik,
  # Lin, Jiang) require `godata(computeIC = TRUE)`, which on this dataset did not complete;
  # Wang is graph-based, reads the GO DAG directly, needs no corpus statistics, and returns
  # in seconds. Passing an IC-based measure here is allowed but will be slow or hang, and
  # `computeIC` is set from `method` accordingly.
  if (!is(enrich_result, "enrichResult")) {
    message("Input is not an enrichResult object; skipping GO reduction.")
    return(enrich_result)
  }
  ontology <- tryCatch(enrich_result@ontology, error = function(e) NULL)
  if (is.null(ontology) || !ontology %in% c("BP", "MF", "CC")) {
    warning("Not a GO enrichment result; skipping GO reduction.")
    return(enrich_result)
  }
  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    warning("GOSemSim not installed; skipping GO reduction and returning the input unchanged.")
    return(enrich_result)
  }
  if (is.null(enrich_result@result) || nrow(enrich_result@result) == 0) {
    warning("Empty enrichment result; skipping GO reduction.")
    return(enrich_result)
  }

  n_original <- nrow(enrich_result@result)
  go_ids <- enrich_result@result$ID

  if (is.null(orgdb)) {
    organism <- tryCatch(enrich_result@organism, error = function(e) "")
    orgdb <- if (grepl("sapiens|human|^hsa$", organism, ignore.case = TRUE)) {
      "org.Hs.eg.db"
    } else if (grepl("musculus|mouse|^mmu$", organism, ignore.case = TRUE)) {
      "org.Mm.eg.db"
    } else {
      warning(sprintf("Cannot map organism '%s' to an OrgDb; skipping GO reduction.", organism))
      return(enrich_result)
    }
  }
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    warning(sprintf("%s not installed; skipping GO reduction.", orgdb))
    return(enrich_result)
  }

  message(sprintf("Reducing GO redundancy: %d terms, %s similarity, threshold %.2f",
                  n_original, method, threshold))

  out <- tryCatch({
    need_ic <- !identical(method, "Wang")
    semdata <- GOSemSim::godata(orgdb, ont = ontology, computeIC = need_ic)

    # Pairwise similarity. O(n^2): fine for the tens-of-terms lists this is meant for,
    # and the reason the step is optional rather than always-on.
    n <- length(go_ids)
    sim <- matrix(0, nrow = n, ncol = n, dimnames = list(go_ids, go_ids))
    for (i in seq_len(n)) {
      sim[i, i] <- 1
      if (i < n) {
        for (j in (i + 1):n) {
          s <- tryCatch(GOSemSim::goSim(go_ids[i], go_ids[j], semData = semdata,
                                        measure = method),
                        error = function(e) 0)
          if (is.na(s)) s <- 0
          sim[i, j] <- sim[j, i] <- s
        }
      }
    }

    # Greedy representative selection: repeatedly take the most significant remaining term
    # and drop everything above `threshold` similarity to it. The representative is chosen
    # by adjusted p-value, so the surviving list is the significant one, not an arbitrary
    # member of each cluster.
    padj <- enrich_result@result$p.adjust
    keep <- integer(0)
    remaining <- seq_len(n)
    while (length(remaining) > 0) {
      best <- remaining[which.min(padj[remaining])]
      keep <- c(keep, best)
      remaining <- remaining[!(sim[best, remaining] > threshold)]
    }
    keep <- sort(keep)

    reduced <- enrich_result
    reduced@result <- enrich_result@result[keep, , drop = FALSE]
    n_reduced <- nrow(reduced@result)
    message(sprintf("GO reduction complete: %d -> %d terms (%.1f%% reduction)",
                    n_original, n_reduced, 100 * (n_original - n_reduced) / n_original))

    attr(reduced, "rrvgo_reduced") <- TRUE
    attr(reduced, "rrvgo_original_n") <- n_original
    attr(reduced, "rrvgo_reduced_n") <- n_reduced
    attr(reduced, "rrvgo_threshold") <- threshold
    attr(reduced, "rrvgo_method") <- method
    reduced
  }, error = function(e) {
    warning(sprintf("GO reduction failed (%s); returning the input unchanged.", e$message))
    enrich_result
  })
  out
}

