test_that("effect_spread = 0 reproduces the default simulation", {
  a <- pe_simulate(seed = 1)
  b <- pe_simulate(effect_spread = 0, seed = 1)
  expect_identical(a, b)
  expect_true(all(a$sim$strength == 1))
  expect_named(a$sim$strength, names(a$sim$targets))
})

test_that("effect_spread varies perturbation strength", {
  x <- pe_simulate(n_single = 40, n_combo = 0, shared_effect = 0,
                   noise_sd = 2, effect_spread = 1, seed = 3)
  st <- x$sim$strength
  expect_gt(stats::sd(log(st)), 0.5)
  eff <- x$sim$effects[, names(st)]
  # effects scale with strength
  expect_gt(stats::cor(sqrt(colSums(eff^2)), st, method = "spearman"), 0.8)
  # and reliability follows strength
  rel <- pe_signal(x, seed = 1)$per_perturbation
  r <- rel$reliability[match(names(st), rel$perturbation)]
  expect_gt(stats::cor(st, r, method = "spearman"), 0.5)
})

test_that("combinations add the scaled single effects", {
  x <- pe_simulate(effect_spread = 1, seed = 2)
  combo <- grep("+", colnames(x$sim$effects), fixed = TRUE, value = TRUE)[1]
  parts <- strsplit(combo, "+", fixed = TRUE)[[1]]
  expect_equal(x$sim$effects[, combo],
               x$sim$effects[, parts[1]] + x$sim$effects[, parts[2]])
})

test_that("effect_spread is validated", {
  expect_error(pe_simulate(effect_spread = -1), "non-negative")
  expect_error(pe_simulate(effect_spread = NA), "non-negative")
})
