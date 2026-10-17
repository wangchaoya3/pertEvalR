test_that("pe_meta_metrics computes DRF and BDS from bounds", {
  b <- data.frame(
    perturbation = rep(c("a", "b", "c", "d"), 2),
    metric = rep(c("pearson_delta", "mse"), each = 4),
    lower = c(0.2, 0.2, 0.4, 0.5, 0.4, 0.4, 0.2, 0.1),
    upper = c(0.8, 0.6, 0.4, 0.2, 0.1, 0.1, 0.2, 0.3)
  )
  mm <- pe_meta_metrics(b)
  pd <- mm[mm$metric == "pearson_delta", ]
  expect_equal(pd$lower, 0.325)
  expect_equal(pd$upper, 0.5)
  expect_equal(pd$ideal, 1)
  expect_equal(pd$drf, (0.5 - 0.325) / (1 - 0.325))
  expect_equal(pd$bds, (1 + 1 + 0.5 + 0) / 4)
  ms <- mm[mm$metric == "mse", ]
  expect_equal(ms$ideal, 0)
  expect_equal(ms$drf, (0.175 - 0.275) / (0 - 0.275))
  expect_equal(ms$bds, (1 + 1 + 0.5 + 0) / 4)
})

test_that("pe_meta_metrics works on pe_bounds output and subsets", {
  x <- pe_simulate(noise_sd = 0.3, seed = 1)
  s <- pe_split(x, seed = 1)
  b <- pe_bounds(x, s, seed = 1)
  mm <- pe_meta_metrics(b)
  expect_setequal(mm$metric, unique(b$metric))
  expect_true(all(mm$n == length(s$test)))
  # with little noise the replicate beats the uninformative baselines
  expect_true(all(mm$bds[mm$metric != "pds"] > 0.9))
  expect_true(all(mm$drf > 0 & mm$drf <= 1.05, na.rm = TRUE))
  sub <- pe_meta_metrics(pe_subset(b, s$test[1:2]))
  expect_true(all(sub$n == 2))
})

test_that("pe_meta_metrics handles missing bounds and validates input", {
  b <- data.frame(perturbation = c("a", "b"), metric = "pearson_delta",
                  lower = NA_real_, upper = c(0.5, 0.6))
  mm <- pe_meta_metrics(b)
  expect_equal(mm$n, 0)
  expect_true(is.na(mm$drf) && is.na(mm$bds))
  expect_error(pe_meta_metrics(list()), "pe_bounds")
})
