test_that("pe_prediction reorders columns to match the split", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  truth <- pe_truth(x, s)
  shuffled <- truth[, rev(colnames(truth)), drop = FALSE]

  p <- pe_prediction(shuffled, s, method = "mine")
  expect_s3_class(p, "pe_prediction")
  expect_identical(colnames(p$pred), s$test)
  expect_equal(p$pred, truth)
  expect_output(print(p), "pe_prediction: mine")
})

test_that("pe_prediction validates its inputs", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  truth <- pe_truth(x, s)

  expect_error(pe_prediction(truth, list()), "pe_split")
  expect_error(pe_prediction(unname(truth), s), "row names")
  expect_error(pe_prediction(truth[, -1, drop = FALSE], s), "held-out")
  expect_error(pe_prediction(truth, s, method = c("a", "b")), "single string")
})
