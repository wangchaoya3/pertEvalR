test_that("unseen_single split is disjoint and excludes the control", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, type = "unseen_single", test_frac = 0.25, seed = 3)

  expect_s3_class(s, "pe_split")
  expect_length(intersect(s$train, s$test), 0)
  expect_false(x$control %in% c(s$train, s$test))
  expect_true(all(s$test %in% pe_perturbations(x, type = "single")))
  expect_length(s$test, 5) # 25% of 20 singles
})

test_that("unseen_single split drops combinations involving held-out genes", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, type = "unseen_single", test_frac = 0.5, seed = 3)
  train_genes <- unlist(strsplit(s$train, "+", fixed = TRUE))
  expect_false(any(train_genes %in% s$test))
})

test_that("unseen_combo split holds out only combinations with known parts", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, type = "unseen_combo", seed = 3)
  singles <- pe_perturbations(x, type = "single")

  expect_true(all(s$test %in% pe_perturbations(x, type = "combo")))
  test_genes <- unlist(strsplit(s$test, "+", fixed = TRUE))
  expect_true(all(test_genes %in% singles))
  expect_true(all(singles %in% s$train))
})

test_that("splits are reproducible and leave the global RNG untouched", {
  x <- pe_simulate(seed = 1)
  set.seed(99)
  before <- .Random.seed
  s1 <- pe_split(x, seed = 7)
  expect_identical(.Random.seed, before)
  s2 <- pe_split(x, seed = 7)
  expect_identical(s1$test, s2$test)
})

test_that("pe_split checks its arguments", {
  x <- pe_simulate(n_combo = 0, seed = 1)
  expect_error(pe_split(x, test_frac = 1), "between 0 and 1")
  expect_error(pe_split(x, type = "unseen_combo"), "No combinatorial")
})

test_that("unseen_single always keeps at least one training single", {
  x <- pe_simulate(n_single = 2, n_combo = 0, seed = 1)
  s <- pe_split(x, test_frac = 0.9, seed = 1)
  expect_length(s$test, 1)
  expect_length(s$train, 1)
})
