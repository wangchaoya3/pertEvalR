all_metrics <- c("mse", "mse_top", "pearson_delta", "pearson_delta_top",
                 "pearson_delta_pert", "pds")

values_of <- function(res, metric) res$value[res$metric == metric]

test_that("pe_evaluate returns one row per perturbation and metric", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  res <- pe_evaluate(baseline_mean(x, s), x)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("method", "perturbation", "metric", "value"))
  expect_identical(nrow(res), length(s$test) * length(all_metrics))
  expect_setequal(unique(res$metric), all_metrics)
  expect_setequal(unique(res$perturbation), s$test)
})

test_that("a perfect prediction gets ideal scores", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  oracle <- pe_prediction(pe_truth(x, s), s, method = "oracle")
  n <- length(s$test)

  for (dist in c("l1", "l2", "cosine")) {
    res <- pe_evaluate(oracle, x, distance = dist)
    expect_equal(values_of(res, "mse"), rep(0, n))
    expect_equal(values_of(res, "mse_top"), rep(0, n))
    expect_equal(values_of(res, "pearson_delta"), rep(1, n))
    expect_equal(values_of(res, "pearson_delta_top"), rep(1, n))
    expect_equal(values_of(res, "pearson_delta_pert"), rep(1, n))
    expect_equal(values_of(res, "pds"), rep(1, n))
  }
})

test_that("mse matches a direct calculation", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_mean(x, s)
  res <- pe_evaluate(p, x, metrics = "mse")
  expected <- unname(colMeans((p$pred - pe_truth(x, s))^2))
  expect_equal(res$value, expected)
})

test_that("constant baselines are undefined on their own reference", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  res_ctrl <- pe_evaluate(baseline_control(x, s), x)
  res_mean <- pe_evaluate(baseline_mean(x, s), x)

  expect_true(all(is.na(values_of(res_ctrl, "pearson_delta"))))
  expect_true(all(is.na(values_of(res_ctrl, "pearson_delta_top"))))
  expect_true(all(is.na(values_of(res_mean, "pearson_delta_pert"))))
  expect_false(anyNA(values_of(res_mean, "pearson_delta")))
})

test_that("constant baselines score at chance on the discrimination score", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, test_frac = 0.4, seed = 2)
  for (p in list(baseline_control(x, s), baseline_mean(x, s))) {
    for (dist in c("l1", "l2", "cosine")) {
      pds <- values_of(pe_evaluate(p, x, metrics = "pds", distance = dist),
                       "pds")
      expect_equal(mean(pds), 0.5)
    }
  }
})

test_that("strong systematic variation flatters the mean baseline", {
  # With a large shift shared by all perturbed cells, the mean baseline gets a
  # high correlation against control but carries no perturbation-specific
  # signal, which the centroid-referenced metric exposes.
  x <- pe_simulate(shared_effect = 2, effect_sd = 0.5, seed = 1)
  s <- pe_split(x, seed = 2)
  res <- pe_evaluate(baseline_mean(x, s), x)
  expect_gt(mean(values_of(res, "pearson_delta")), 0.8)
  expect_true(all(is.na(values_of(res, "pearson_delta_pert"))))
})

test_that("pe_evaluate accepts a list of predictions", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  res <- pe_evaluate(list(baseline_control(x, s), baseline_mean(x, s)), x,
                     metrics = c("mse", "pds"))
  expect_setequal(unique(res$method), c("control", "mean"))
  expect_identical(nrow(res), 2L * 2L * length(s$test))
})

test_that("pe_evaluate checks its arguments", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_mean(x, s)

  expect_error(pe_evaluate(p$pred, x), "pe_prediction")
  expect_error(pe_evaluate(p, x, metrics = "nonsense"))
  expect_error(pe_evaluate(p, x, top_n = 0), "positive")

  p_missing <- p
  p_missing$pred <- p_missing$pred[-1, , drop = FALSE]
  expect_error(pe_evaluate(p_missing, x), "lacks")
})

test_that("top_n larger than the number of genes is capped", {
  x <- pe_simulate(n_genes = 30, n_single = 5, n_combo = 0,
                   n_affected = 5, seed = 1)
  s <- pe_split(x, seed = 2)
  p <- baseline_mean(x, s)
  res_big <- pe_evaluate(p, x, metrics = "mse_top", top_n = 1000)
  res_all <- pe_evaluate(p, x, metrics = "mse")
  expect_equal(res_big$value, res_all$value)
})

test_that("pds is NA with a single held-out perturbation", {
  x <- pe_simulate(n_single = 3, n_combo = 0, seed = 1)
  s <- pe_split(x, test_frac = 0.2, seed = 2)
  expect_length(s$test, 1)
  res <- pe_evaluate(baseline_mean(x, s), x, metrics = "pds")
  expect_true(is.na(res$value))
})

test_that("pe_summary averages over perturbations", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 2)
  res <- pe_evaluate(list(baseline_control(x, s), baseline_mean(x, s)), x)
  sm <- pe_summary(res)

  expect_named(sm, c("method", "metric", "mean", "sd", "n", "n_na"))
  expect_identical(nrow(sm), 2L * length(all_metrics))
  expect_true(all(sm$n + sm$n_na == length(s$test)))

  row <- sm$method == "mean" & sm$metric == "mse"
  expect_equal(sm$mean[row], mean(res$value[res$method == "mean" &
                                              res$metric == "mse"]))
  na_row <- sm$method == "control" & sm$metric == "pearson_delta"
  expect_true(is.na(sm$mean[na_row]))
  expect_identical(sm$n[na_row], 0L)

  expect_error(pe_summary(data.frame(a = 1)), "pe_evaluate")
})
