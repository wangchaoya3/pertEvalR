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

# baseline_additive ---------------------------------------------------------

test_that("baseline_additive adds the single-perturbation changes", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, type = "unseen_combo", seed = 2)
  p <- baseline_additive(x, s)

  expect_identical(p$method, "additive")
  expect_identical(colnames(p$pred), s$test)

  combo <- s$test[1]
  parts <- strsplit(combo, "+", fixed = TRUE)[[1]]
  pb <- pe_pseudobulk(x, c("control", parts))
  expected <- pb[, "control"] +
    (pb[, parts[1]] - pb[, "control"]) +
    (pb[, parts[2]] - pb[, "control"])
  expect_equal(p$pred[, combo], expected)
})

test_that("baseline_additive is close to the truth for additive simulations", {
  # With no shared shift the simulated combinations are exactly additive.
  # (A shared shift is counted once per component by the additive baseline.)
  x <- pe_simulate(noise_sd = 0.1, shared_effect = 0, seed = 1)
  s <- pe_split(x, type = "unseen_combo", seed = 2)
  truth <- pe_truth(x, s)
  mse_add <- mean((baseline_additive(x, s)$pred - truth)^2)
  mse_mean <- mean((baseline_mean(x, s)$pred - truth)^2)
  expect_lt(mse_add, mse_mean)
})

test_that("baseline_additive rejects splits without combinations", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, type = "unseen_single", seed = 2)
  expect_error(baseline_additive(x, s), "combinatorial perturbations only")
})

# baseline_linear -----------------------------------------------------------

test_that("baseline_linear returns finite predictions of the right shape", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_linear(x, s, targets = x$sim$targets)

  expect_identical(p$method, "linear")
  expect_identical(dim(p$pred), dim(pe_truth(x, s)))
  expect_identical(colnames(p$pred), s$test)
  expect_true(all(is.finite(p$pred)))
})

test_that("baseline_linear tends to the training mean as ridge grows", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_linear(x, s, ridge = 1e12, targets = x$sim$targets)

  train_singles <- s$train[!grepl("+", s$train, fixed = TRUE)]
  b <- rowMeans(pe_pseudobulk(x, train_singles))
  for (j in seq_len(ncol(p$pred))) {
    expect_equal(p$pred[, j], b)
  }
})

test_that("baseline_linear caps n_comp at the number of training perturbations", {
  x <- pe_simulate(n_single = 4, n_combo = 0, seed = 1)
  s <- pe_split(x, test_frac = 0.25, seed = 2)
  expect_no_error(baseline_linear(x, s, n_comp = 50, targets = x$sim$targets))
})

test_that("baseline_linear checks targets and split type", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  # Simulated labels (P1, P2, ...) are not gene names.
  expect_error(baseline_linear(x, s), "Target genes not found")
  expect_error(baseline_linear(x, s, targets = c(P1 = "g1")), "no entry")
  expect_error(baseline_linear(x, s, targets = unname(x$sim$targets)),
               "named character vector")

  s_combo <- pe_split(x, type = "unseen_combo", seed = 2)
  expect_error(baseline_linear(x, s_combo, targets = x$sim$targets),
               "single perturbations only")
})
