test_that("pe_subset keeps rows and settings", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 1)
  res <- pe_evaluate(list(baseline_mean(x, s), baseline_control(x, s)), x,
                     top_n = 10)
  keep <- s$test[1:2]
  sub <- pe_subset(res, keep)
  expect_setequal(unique(sub$perturbation), keep)
  expect_equal(nrow(sub), sum(res$perturbation %in% keep))
  expect_identical(attr(sub, "settings"), attr(res, "settings"))
  expect_identical(rownames(sub), as.character(seq_len(nrow(sub))))
})

test_that("pe_subset output calibrates against subset bounds", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 1)
  res <- pe_evaluate(baseline_mean(x, s), x, top_n = 10)
  b <- pe_bounds(x, s, top_n = 10, seed = 1)
  keep <- s$test[1:2]
  cal <- pe_calibrate(pe_subset(res, keep), pe_subset(b, keep))
  # metrics the mean baseline cannot score (NA) have n = 0
  expect_true(all(cal$n[!is.na(cal$score)] == 2))
  expect_true(any(!is.na(cal$score)))
  # settings survive, so a mismatch is still caught
  b20 <- pe_bounds(x, s, top_n = 20, seed = 1)
  expect_error(pe_calibrate(pe_subset(res, keep), pe_subset(b20, keep)),
               "different")
  # subsetting only one side is caught too
  expect_error(pe_calibrate(pe_subset(res, keep), b), "different perturbations")
})

test_that("pe_subset accepts a pe_signal object", {
  x <- pe_simulate(seed = 1)
  s <- pe_split(x, seed = 1)
  res <- pe_evaluate(baseline_mean(x, s), x)
  sig <- pe_signal(x, seed = 1)
  sig$per_perturbation$reliability <- 0
  sig$per_perturbation$reliability[sig$per_perturbation$perturbation ==
                                     s$test[1]] <- 0.9
  sub <- pe_subset(res, sig, min_reliability = 0.5)
  expect_identical(unique(sub$perturbation), s$test[1])
  expect_error(pe_subset(res, sig, min_reliability = 0.95), "None")
  expect_error(pe_subset(res, sig, min_reliability = NA), "single number")
})

test_that("pe_subset validates its input", {
  expect_error(pe_subset(list(), "a"), "data frame")
  df <- data.frame(perturbation = c("a", "b"), value = 1:2)
  expect_error(pe_subset(df, 1), "character")
  expect_error(pe_subset(df, "z"), "None")
})
