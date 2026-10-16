test_that("baseline_control predicts the control mean for every test perturbation", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_control(x, s)

  expect_s3_class(p, "pe_prediction")
  expect_identical(p$method, "control")
  expect_identical(colnames(p$pred), s$test)
  expect_identical(rownames(p$pred), rownames(x$expr))

  ctrl_mean <- rowMeans(x$expr[, x$perturbation == "control"])
  for (j in seq_len(ncol(p$pred))) {
    expect_equal(p$pred[, j], ctrl_mean)
  }
})

test_that("baseline_mean averages training pseudobulks with equal weights", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_mean(x, s)

  per_pert <- vapply(s$train,
                     function(l) rowMeans(x$expr[, x$perturbation == l]),
                     numeric(nrow(x$expr)))
  expect_equal(p$pred[, 1], rowMeans(per_pert))
  expect_identical(dim(p$pred), dim(pe_truth(x, s)))
})

test_that("baseline_mean never uses held-out perturbations", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  # Corrupt the held-out cells; the prediction must not change.
  y <- x
  y$expr[, y$perturbation %in% s$test] <- 1e6
  expect_equal(baseline_mean(x, s)$pred, baseline_mean(y, s)$pred)
})

test_that("with strong systematic variation the mean baseline beats no-change", {
  x <- pe_simulate(shared_effect = 2, effect_sd = 0.5, seed = 1)
  s <- pe_split(x, seed = 2)
  truth <- pe_truth(x, s)
  mse_ctrl <- mean((baseline_control(x, s)$pred - truth)^2)
  mse_mean <- mean((baseline_mean(x, s)$pred - truth)^2)
  expect_lt(mse_mean, mse_ctrl)
})

test_that("baselines reject a split made for a different dataset", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  s$test <- c(s$test, "not-a-perturbation")
  expect_error(baseline_control(x, s), "not present")
})
