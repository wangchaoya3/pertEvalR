all_metrics <- c("mse", "mse_top", "pearson_delta", "pearson_delta_top",
                 "pearson_delta_pert", "pds")

calib_fixture <- function() {
  x <- pe_simulate(n_single = 30, n_combo = 0, cells_per_pert = 40,
                   noise_sd = 0.3, seed = 1)
  s <- pe_split(x, test_frac = 0.3, seed = 2)
  list(x = x, s = s)
}

# pe_bounds ------------------------------------------------------------------

test_that("pe_bounds returns one row per perturbation and metric", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, seed = 1)

  expect_named(b, c("perturbation", "metric", "lower", "lower_method",
                    "upper", "upper_raw"))
  expect_identical(nrow(b), length(f$s$test) * length(all_metrics))
  expect_setequal(unique(b$perturbation), f$s$test)
  expect_identical(attr(b, "settings"), list(top_n = 20L, distance = "l1"))
})

test_that("the lower bound is the better uninformative baseline", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, seed = 1)
  base <- pe_evaluate(list(baseline_control(f$x, f$s),
                           baseline_mean(f$x, f$s)), f$x)

  # pearson_delta: undefined for no-change, so the mean baseline sets it
  pd <- b[b$metric == "pearson_delta", ]
  expect_true(all(pd$lower_method == "mean"))
  mean_pd <- base[base$method == "mean" & base$metric == "pearson_delta", ]
  expect_equal(pd$lower, mean_pd$value[match(pd$perturbation,
                                             mean_pd$perturbation)])

  # pearson_delta_pert: undefined for the mean, so no-change sets it
  pp <- b[b$metric == "pearson_delta_pert", ]
  expect_true(all(pp$lower_method == "control"))

  # mse: the baseline with the smaller average error, used for every
  # perturbation (one choice per metric, not per perturbation)
  ms <- b[b$metric == "mse", ]
  expect_length(unique(ms$lower_method), 1)
  ctrl <- base[base$method == "control" & base$metric == "mse", ]
  mn <- base[base$method == "mean" & base$metric == "mse", ]
  best <- if (mean(ctrl$value) <= mean(mn$value)) ctrl else mn
  expect_identical(ms$lower_method[1], best$method[1])
  expect_equal(ms$lower, best$value[match(ms$perturbation, best$perturbation)])
})

test_that("the lower bound is not inflated by per-perturbation selection", {
  # Both baselines make one prediction for every perturbation, so each has a
  # mean discrimination score of exactly 0.5; so must the lower bound.
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, metrics = "pds", seed = 1)
  expect_equal(mean(b$lower), 0.5)
})

test_that("the upper bound is corrected from half to full sample size", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, seed = 1)

  for (m in c("mse", "mse_top")) {
    d <- b[b$metric == m, ]
    expect_equal(d$upper, d$upper_raw / 4)
  }
  for (m in c("pearson_delta", "pearson_delta_top", "pearson_delta_pert")) {
    d <- b[b$metric == m, ]
    r <- d$upper_raw
    expected <- ifelse(r > 0, sqrt(2 * r / (1 + r)), r)
    expect_equal(d$upper, expected)
    expect_true(all(d$upper >= d$upper_raw))
    expect_true(all(d$upper <= 1))
  }
  d <- b[b$metric == "pds", ]
  expect_equal(d$upper, d$upper_raw)
})

test_that("the replicate upper bound beats the uninformative lower bound", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, seed = 1)
  agg <- function(m, col) mean(b[[col]][b$metric == m])

  for (m in c("pearson_delta", "pearson_delta_top", "pearson_delta_pert",
              "pds")) {
    expect_gt(agg(m, "upper"), agg(m, "lower"))
  }
  for (m in c("mse", "mse_top")) {
    expect_lt(agg(m, "upper"), agg(m, "lower"))
  }
  expect_gt(agg("pearson_delta", "upper"), 0.8)
})

test_that("pe_bounds is reproducible and leaves the global RNG untouched", {
  f <- calib_fixture()
  set.seed(11)
  before <- .Random.seed
  b1 <- pe_bounds(f$x, f$s, metrics = "mse", seed = 5)
  expect_identical(.Random.seed, before)
  b2 <- pe_bounds(f$x, f$s, metrics = "mse", seed = 5)
  expect_identical(b1, b2)
})

test_that("pe_bounds needs at least two cells per held-out perturbation", {
  expr <- matrix(rnorm(30), nrow = 3,
                 dimnames = list(paste0("g", 1:3), NULL))
  pert <- c(rep("control", 4), "A", "B", "B", "C", "C", "C")
  x <- pe_data(expr, pert)
  s <- structure(list(train = c("B", "C"), test = "A", control = "control",
                      type = "unseen_single", test_frac = 0.3, seed = NULL),
                 class = "pe_split")
  expect_error(pe_bounds(x, s), "fewer than two cells")
})

# pe_calibrate --------------------------------------------------------------

test_that("pe_calibrate places scores between the bounds", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, seed = 1)
  oracle <- pe_prediction(pe_truth(f$x, f$s), f$s, method = "oracle")
  res <- pe_evaluate(list(oracle, baseline_mean(f$x, f$s)), f$x)
  cal <- pe_calibrate(res, b)

  expect_named(cal, c("method", "metric", "n", "score", "lower", "upper",
                      "calibrated", "calibrated_median"))
  expect_identical(nrow(cal), 2L * length(all_metrics))

  # The oracle is compared with the very data it copies, so it reaches or
  # exceeds the replicate bound on every metric.
  or <- cal[cal$method == "oracle", ]
  expect_true(all(or$calibrated >= 1 - 1e-8))

  # An uninformative baseline cannot beat the better uninformative baseline,
  # and scores exactly 0 on the metrics where it is that baseline.
  mn <- cal[cal$method == "mean" & !is.na(cal$calibrated), ]
  expect_true(nrow(mn) > 0)
  expect_true(all(mn$calibrated <= 1e-8))
  mean_is_lower <- unique(b$metric[b$lower_method %in% "mean"])
  expect_equal(mn$calibrated[mn$metric %in% mean_is_lower],
               rep(0, sum(mn$metric %in% mean_is_lower)))
})

test_that("calibrated equals the formula on the averaged bounds", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, metrics = "mse", seed = 1)
  res <- pe_evaluate(baseline_linear(f$x, f$s, targets = f$x$sim$targets),
                     f$x, metrics = "mse")
  cal <- pe_calibrate(res, b)
  expect_equal(cal$calibrated, (cal$score - cal$lower) / (cal$upper - cal$lower))
  expect_equal(cal$score, mean(res$value))
})

test_that("metrics without headroom are not calibrated", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, metrics = "mse", seed = 1)
  b$upper <- b$lower           # no room between the bounds
  res <- pe_evaluate(baseline_mean(f$x, f$s), f$x, metrics = "mse")
  cal <- pe_calibrate(res, b)
  expect_true(is.na(cal$calibrated))
  expect_true(is.na(cal$calibrated_median))
})

test_that("undefined scores give an empty calibration row", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, metrics = "pearson_delta", seed = 1)
  res <- pe_evaluate(baseline_control(f$x, f$s), f$x,
                     metrics = "pearson_delta")
  cal <- pe_calibrate(res, b)
  expect_identical(cal$n, 0L)
  expect_true(is.na(cal$calibrated))
})

test_that("pe_calibrate checks that results and bounds match", {
  f <- calib_fixture()
  b <- pe_bounds(f$x, f$s, metrics = "mse", seed = 1)
  res <- pe_evaluate(baseline_mean(f$x, f$s), f$x, metrics = "mse")

  expect_error(pe_calibrate(data.frame(a = 1), b), "pe_evaluate")
  expect_error(pe_calibrate(res, data.frame(a = 1)), "pe_bounds")

  other <- pe_split(f$x, test_frac = 0.3, seed = 99)
  res_other <- pe_evaluate(baseline_mean(f$x, other), f$x, metrics = "mse")
  expect_error(pe_calibrate(res_other, b), "different perturbations")

  res_pds <- pe_evaluate(baseline_mean(f$x, f$s), f$x, metrics = "pds")
  expect_error(pe_calibrate(res_pds, b), "lacks metric")

  res_top <- pe_evaluate(baseline_mean(f$x, f$s), f$x, metrics = "mse",
                         top_n = 50)
  expect_error(pe_calibrate(res_top, b), "different")
})
