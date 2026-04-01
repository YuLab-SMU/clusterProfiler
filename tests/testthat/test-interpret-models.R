library(testthat)

make_interpret_input <- function() {
  data.frame(
    ID = c("GO:0006915", "GO:0008284"),
    Description = c("apoptotic process", "positive regulation of proliferation"),
    GeneRatio = c("10/100", "20/100"),
    p.adjust = c(0.01, 0.02),
    geneID = c("TP53/BAX", "MYC/CCND1/CDK4"),
    stringsAsFactors = FALSE
  )
}

test_that("normalize_interpret_model supports NULL, bare IDs, and model objects", {
  normalize_interpret_model <- getFromNamespace(".normalize_interpret_model", "clusterProfiler")

  expect_null(normalize_interpret_model(NULL))

  expect_warning(
    expect_equal(
      normalize_interpret_model("deepseek-chat"),
      "deepseek:deepseek-chat"
    ),
    "Inferred model ID"
  )

  model_obj <- structure(list(provider = "openai", model_id = "gpt-4o-mini"), class = "LanguageModelV1")
  expect_identical(normalize_interpret_model(model_obj), model_obj)
})

test_that("interpret passes NULL model through to aisdk default resolution", {
  captured_model <- "unset"

  local_mocked_bindings(
    .call_generate_object = function(model, ...) {
      captured_model <<- model
      structure(list(
        overview = "Synthetic interpretation",
        confidence = "High"
      ), class = c("interpretation", "list"))
    },
    .package = "clusterProfiler"
  )

  res <- clusterProfiler::interpret(make_interpret_input(), model = NULL)

  expect_null(captured_model)
  expect_s3_class(res, "interpretation")
  expect_equal(res$overview, "Synthetic interpretation")
})

test_that("interpret normalizes bare model IDs before dispatch", {
  captured_model <- NULL

  local_mocked_bindings(
    .call_generate_object = function(model, ...) {
      captured_model <<- model
      structure(list(
        overview = "Synthetic interpretation",
        confidence = "High"
      ), class = c("interpretation", "list"))
    },
    .package = "clusterProfiler"
  )

  expect_warning(
    clusterProfiler::interpret(make_interpret_input(), model = "deepseek-chat"),
    "Inferred model ID"
  )

  expect_equal(captured_model, "deepseek:deepseek-chat")
})

test_that("interpret wrappers default to aisdk global model via model = NULL", {
  expect_null(formals(clusterProfiler::interpret)$model)
  expect_null(formals(clusterProfiler::interpret_agent)$model)
  expect_null(formals(clusterProfiler::interpret_hierarchical)$model)
})
