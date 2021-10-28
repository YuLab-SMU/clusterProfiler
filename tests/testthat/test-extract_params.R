context("Test extract_params")

test_function = function(x, y, z = 'z'){
  this_call = match.call(expand.dots = TRUE)
  return(extract_params(this_call))
}

test_that("Can extract position parameters", {
  params = test_function(1, 2, 3)
  expect_equal(params, list(x = "1", y = "2", z = "3"))
  # note these are coerced to character, might not be desirable..
})

test_that("Can extract keyword parameters", {
  params = test_function(x = 1, z = 3, y = 2)
  expect_equal(params, list(x = "1", y = "2", z = "3"))
})

test_that("Can extract default parameters", {
  params = test_function(x = 1, y = 2)
  expect_equal(params, list(x = "1", y = "2", z = "z"))
})
