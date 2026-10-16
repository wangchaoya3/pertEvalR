compare_fixture <- function() {
  x <- pe_simulate(n_single = 40, n_combo = 0, seed = 1)
  s <- pe_split(x, test_frac = 0.5, seed = 2)   # 20 held-out perturbations
  oracle <- pe_prediction(pe_truth(x, s), s, method = "oracle")
  list(x = x, s = s,
       preds = list(oracle, baseline_mean(x, s), baseline_control(x, s)))
}

test_that("pe_compare returns one row per metric and pair", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds, f$x, metrics = c("mse", "pds"))
  cmp <- pe_compare(res, n_boot = 200, seed = 1)

  expect_named(cmp, c("metric", "method_a", "method_b", "n", "mean_a",
                      "mean_b", "diff", "ci_low", "ci_high", "p_value",
                      "p_adj", "favours"))
  expect_identical(nrow(cmp), 2L * 3L)   # 2 metrics x choose(3, 2) pairs
  expect_true(all(cmp$n == 20))
})

test_that("reference mode compares every method with the reference", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds, f$x, metrics = "mse")
  cmp <- pe_compare(res, reference = "mean", n_boot = 200, seed = 1)

  expect_identical(nrow(cmp), 2L)
  expect_true(all(cmp$method_b == "mean"))
  expect_setequal(cmp$method_a, c("oracle", "control"))
})

test_that("a clearly better method is favoured, respecting metric direction", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds[1:2], f$x, metrics = c("mse", "pds"))
  cmp <- pe_compare(res, n_boot = 500, seed = 1)

  mse <- cmp[cmp$metric == "mse", ]
  expect_lt(mse$diff, 0)            # oracle has lower error
  expect_lt(mse$ci_high, 0)
  expect_identical(mse$favours, "oracle")
  expect_lt(mse$p_value, 0.001)

  pds <- cmp[cmp$metric == "pds", ]
  expect_gt(pds$diff, 0)            # oracle discriminates better
  expect_identical(pds$favours, "oracle")
})

test_that("diff, means and interval agree with a direct calculation", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds[1:2], f$x, metrics = "mse")
  cmp <- pe_compare(res, n_boot = 500, seed = 1)

  a <- res$value[res$method == "oracle"]
  b <- res$value[res$method == "mean"]
  expect_equal(cmp$mean_a, mean(a))
  expect_equal(cmp$mean_b, mean(b))
  expect_equal(cmp$diff, mean(a - b))
  expect_lte(cmp$ci_low, cmp$diff)
  expect_gte(cmp$ci_high, cmp$diff)
})

test_that("identical methods are not distinguished", {
  f <- compare_fixture()
  p1 <- baseline_mean(f$x, f$s)
  p2 <- p1
  p2$method <- "mean_copy"
  res <- pe_evaluate(list(p1, p2), f$x, metrics = "mse")
  cmp <- pe_compare(res, n_boot = 200, seed = 1)

  expect_equal(cmp$diff, 0)
  expect_equal(c(cmp$ci_low, cmp$ci_high), c(0, 0))
  expect_equal(cmp$p_value, 1)
  expect_identical(cmp$favours, "neither")
})

test_that("undefined metric values are dropped from the comparison", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds, f$x, metrics = "pearson_delta")
  cmp <- pe_compare(res, reference = "mean", n_boot = 200, seed = 1)

  # pearson_delta is NA for every perturbation of the no-change baseline
  ctrl <- cmp[cmp$method_a == "control", ]
  expect_identical(ctrl$n, 0L)
  expect_true(is.na(ctrl$diff))
  expect_true(is.na(ctrl$favours))
})

test_that("p values are adjusted across rows", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds, f$x, metrics = c("mse", "pds"))
  cmp <- pe_compare(res, n_boot = 200, seed = 1)
  ok <- !is.na(cmp$p_value)
  expect_true(all(cmp$p_adj[ok] >= cmp$p_value[ok]))

  none <- pe_compare(res, n_boot = 200, seed = 1, p_adjust = "none")
  expect_equal(none$p_adj, none$p_value)
})

test_that("intervals are reproducible and leave the global RNG untouched", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds[1:2], f$x, metrics = "mse")
  set.seed(42)
  before <- .Random.seed
  c1 <- pe_compare(res, n_boot = 200, seed = 7)
  expect_identical(.Random.seed, before)
  c2 <- pe_compare(res, n_boot = 200, seed = 7)
  expect_identical(c1, c2)
})

test_that("pe_compare checks its arguments", {
  f <- compare_fixture()
  res <- pe_evaluate(f$preds[1:2], f$x, metrics = "mse")

  expect_error(pe_compare(data.frame(a = 1)), "pe_evaluate")
  expect_error(pe_compare(rbind(res, res)), "repeated")
  expect_error(pe_compare(res, metrics = "pds"), "not found")
  expect_error(pe_compare(res[res$method == "mean", ]), "two methods")
  expect_error(pe_compare(res, reference = "nope"), "reference")
  expect_error(pe_compare(res, n_boot = 0), "positive")
  expect_error(pe_compare(res, conf = 1), "between 0 and 1")
  expect_error(pe_compare(res, p_adjust = "magic"))
})
