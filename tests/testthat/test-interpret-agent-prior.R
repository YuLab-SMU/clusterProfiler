library(testthat)

agent_input <- function() {
  data.frame(
    ID = c("GO:0042110", "GO:0002250"),
    Description = c("T cell activation", "adaptive immune response"),
    GeneRatio = c("10/100", "8/100"),
    p.adjust = c(0.001, 0.01),
    geneID = c("CD3E/CD4/IL7R/CCR7", "CD3E/CD4/LCK"),
    stringsAsFactors = FALSE
  )
}

# Replace the model call: record every prompt and schema, answer each stage with a minimal object.
mock_agent_calls <- function(env) {
  local_mocked_bindings(
    .gen_object_stamped = function(reasoning, collector, stage, model, prompt, schema, schema_name,
                                   system, temperature, max_tokens) {
      env$calls[[stage]] <- list(prompt = prompt, schema = schema, system = system)
      obj <- switch(stage,
        cleaner = list(kept_pathways = list("T cell activation"), discarded_pathways = list(), reasoning = "r"),
        detective = list(key_drivers = list("TCF7"), functional_modules = list("TCR"), refined_network = list(),
                         network_evidence = "e"),
        synthesizer = list(overview = "CD4 T cells", conclusion = "naive CD4 T cell",
                           prior_verdict = "Refined", prior_rationale = "CCR7 and IL7R"),
        list()
      )
      list(object = obj, finish_reason = "stop", raw_text = "{}")
    },
    .package = "clusterProfiler", .env = parent.frame()
  )
}

test_that("interpret() no longer takes task or prior, and sends every reduced term by default", {
  f <- formals(clusterProfiler::interpret)
  expect_false("task" %in% names(f))
  expect_false("prior" %in% names(f))
  expect_null(f$n_pathways)
  expect_false("task" %in% names(formals(clusterProfiler::interpret_hierarchical)))
  expect_null(formals(clusterProfiler::interpret_hierarchical)$n_pathways)
  expect_true("prior" %in% names(formals(clusterProfiler::interpret_agent)))
  expect_null(formals(clusterProfiler::interpret_agent)$n_pathways)
  expect_error(clusterProfiler::interpret(agent_input(), task = "cell_type"), "unused argument")
})

test_that("interpret_agent() passes prior to the Detective and the Synthesizer as a claim to check", {
  env <- new.env(); env$calls <- list()
  mock_agent_calls(env)
  res <- clusterProfiler::interpret_agent(
    agent_input(), context = "PBMC; which cell type is this cluster?",
    prior = "SingleR (HPCA): T_cell:CD4+_central_memory", model = "openai:gpt-4o-mini"
  )
  for (st in c("detective", "synthesizer")) {
    expect_match(env$calls[[st]]$prompt, "Results of other tools", fixed = TRUE)
    expect_match(env$calls[[st]]$prompt, "T_cell:CD4+_central_memory", fixed = TRUE)
  }
  expect_false(grepl("Results of other tools", env$calls$cleaner$prompt, fixed = TRUE))
  expect_true(all(c("conclusion", "prior_verdict", "prior_rationale") %in% names(env$calls$synthesizer$schema$properties)))
  expect_equal(res$conclusion, "naive CD4 T cell")
  expect_equal(res$prior_verdict, "Refined")
  expect_true(isTRUE(clusterProfiler::provenance(res)$prior_supplied))
})

test_that("without prior the Synthesizer still answers the question but gives no verdict", {
  env <- new.env(); env$calls <- list()
  mock_agent_calls(env)
  res <- clusterProfiler::interpret_agent(agent_input(), context = "which cell type?", model = "openai:gpt-4o-mini")
  props <- names(env$calls$synthesizer$schema$properties)
  expect_true("conclusion" %in% props)
  expect_false("prior_verdict" %in% props)
  expect_false(grepl("Results of other tools", env$calls$synthesizer$prompt, fixed = TRUE))
  expect_false(isTRUE(clusterProfiler::provenance(res)$prior_supplied))
})

test_that("a named prior is matched to its cluster", {
  rp <- getFromNamespace(".resolve_prior", "clusterProfiler")
  expect_equal(rp(c(C1 = "B cell", C2 = "NK cell"), "C2"), "NK cell")
  expect_null(rp(c(C1 = "B cell"), "C3"))
  expect_equal(rp("T cell", "anything"), "T cell")
})
