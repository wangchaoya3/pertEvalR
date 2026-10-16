test_that("pe_signal returns one row per perturbation with cell counts", {
  x <- pe_simulate(seed = 1)
  sig <- pe_signal(x, seed = 1)
  expect_s3_class(sig, "pe_signal")
  per <- sig$per_perturbation
  expect_named(per, c("perturbation", "n_cells", "effect_size",
                      "cosine_to_average", "effective_genes", "reliability"))
  expect_setequal(per$perturbation, pe_perturbations(x))
  expect_true(all(per$n_cells == 30))
  expect_equal(sig$summary$n_perturbations, length(pe_perturbations(x)))
  expect_true(all(per$effect_size > 0))
  expect_true(all(abs(per$cosine_to_average) <= 1 + 1e-12))
})

test_that("strong shared shift gives high cosine and systematic fraction", {
  strong <- pe_signal(pe_simulate(shared_effect = 2, seed = 1), seed = 1)
  weak <- pe_signal(pe_simulate(shared_effect = 0, seed = 1), seed = 1)
  expect_gt(strong$summary$mean_cosine, 0.9)
  expect_gt(strong$summary$systematic_fraction, 0.8)
  expect_lt(weak$summary$mean_cosine, 0.5)
  expect_lt(weak$summary$systematic_fraction, 0.3)
  sf <- c(strong$summary$systematic_fraction, weak$summary$systematic_fraction)
  expect_true(all(sf >= 0 & sf <= 1))
})

test_that("effective_genes reflects sparse effects", {
  x <- pe_simulate(shared_effect = 0, n_affected = 20, seed = 2)
  eg <- pe_signal(x, seed = 1)$per_perturbation$effective_genes
  expect_true(all(eg > 2 & eg < 60))
})

test_that("effective_genes equals the number of equal non-zero changes", {
  genes <- paste0("g", 1:10)
  expr <- matrix(0, 10, 6, dimnames = list(genes, paste0("c", 1:6)))
  expr[1:4, 3:4] <- 1   # A changes four genes by the same amount
  expr[5, 5:6] <- 3     # B changes one gene
  x <- pe_data(expr, c("control", "control", "A", "A", "B", "B"))
  per <- pe_signal(x, seed = 1)$per_perturbation
  expect_equal(per$effective_genes[per$perturbation == "A"], 4)
  expect_equal(per$effective_genes[per$perturbation == "B"], 1)
})

test_that("reliability is high for low noise and NA for single cells", {
  x <- pe_simulate(noise_sd = 0.05, seed = 3)
  rel <- pe_signal(x, seed = 1)$per_perturbation$reliability
  expect_true(all(rel > 0.9))

  set.seed(10)
  expr <- matrix(rnorm(20 * 9), 20, 9,
                 dimnames = list(paste0("g", 1:20), paste0("c", 1:9)))
  y <- pe_data(expr, c(rep("control", 4), "A", rep("B", 4)))
  per <- pe_signal(y, seed = 1)$per_perturbation
  expect_true(is.na(per$reliability[per$perturbation == "A"]))
  expect_false(is.na(per$reliability[per$perturbation == "B"]))
  expect_equal(per$n_cells[per$perturbation == "A"], 1)
})

test_that("pe_signal is reproducible and leaves the global RNG alone", {
  x <- pe_simulate(seed = 1)
  set.seed(99)
  before <- .Random.seed
  a <- pe_signal(x, seed = 5)
  expect_identical(.Random.seed, before)
  b <- pe_signal(x, seed = 5)
  expect_identical(a, b)
})

test_that("a subset of perturbations keeps the dataset-wide average", {
  x <- pe_simulate(seed = 1)
  full <- pe_signal(x, seed = 1)$per_perturbation
  sub <- pe_signal(x, perturbations = c("P1", "P2"), seed = 1)
  expect_equal(sub$per_perturbation$perturbation, c("P1", "P2"))
  expect_equal(sub$per_perturbation$cosine_to_average,
               full$cosine_to_average[match(c("P1", "P2"), full$perturbation)])
  expect_equal(sub$summary$n_perturbations, 2)
})

test_that("pe_signal validates its input", {
  x <- pe_simulate(seed = 1)
  expect_error(pe_signal(list()), "pe_data")
  expect_error(pe_signal(x, perturbations = "nope"), "Unknown")
  expr <- matrix(1, 3, 4, dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  one <- pe_data(expr, c("control", "control", "A", "A"))
  expect_error(pe_signal(one), "at least two")
})

test_that("print.pe_signal reports the systematic fraction", {
  sig <- pe_signal(pe_simulate(seed = 1), seed = 1)
  expect_output(print(sig), "systematic fraction")
})

test_that("pe_bounds still works with the shared halving helper", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 1)
  b <- pe_bounds(x, s, metrics = "pearson_delta", seed = 1)
  expect_true(all(is.finite(b$upper)))
  expect_true(all(b$upper >= b$upper_raw))
})
