#' Compare methods statistically
#'
#' Paired comparison of methods evaluated on the same held-out perturbations.
#' For each metric and pair of methods, the per-perturbation difference
#' `a - b` is summarised by its mean, a percentile bootstrap confidence
#' interval (resampling perturbations), and a paired Wilcoxon signed-rank
#' test. Existing evaluation tools mostly report point estimates; small
#' differences between models are often within resampling noise.
#'
#' Perturbations where either method's value is `NA` are dropped from that
#' comparison (for example, `pearson_delta` is undefined for
#' [baseline_control()]). The `favours` column names the better method when
#' the confidence interval excludes zero, taking into account that lower is
#' better for `mse` and `mse_top` and higher is better for the other metrics;
#' otherwise it is `"neither"`.
#'
#' Inference treats the held-out perturbations as independent samples. That
#' is approximate for `"pds"`, whose value for one perturbation depends on the
#' others held out with it, and it says nothing about variation between data
#' splits; repeat the evaluation over several splits to capture that.
#'
#' @param results A data frame returned by [pe_evaluate()], covering at least
#'   two methods evaluated on the same split. Each combination of method,
#'   perturbation and metric must appear once.
#' @param metrics Metrics to compare. Defaults to all metrics in `results`.
#' @param reference Optional method name. If given, every other method is
#'   compared with it (`method_b` is the reference); otherwise all pairs of
#'   methods are compared.
#' @param n_boot Number of bootstrap resamples.
#' @param conf Confidence level of the bootstrap interval.
#' @param p_adjust Multiple-testing correction applied across all rows, passed
#'   to [stats::p.adjust()] (for example `"holm"`, `"BH"` or `"none"`).
#' @param seed Optional integer seed for reproducible intervals. The global
#'   random number generator state is restored afterwards.
#'
#' @return A data frame with one row per metric and pair of methods and
#'   columns `metric`, `method_a`, `method_b`, `n` (perturbations compared),
#'   `mean_a`, `mean_b`, `diff` (mean of `a - b`), `ci_low`, `ci_high`,
#'   `p_value`, `p_adj` and `favours`.
#' @export
#' @examples
#' x <- pe_simulate(n_single = 40, n_combo = 0, seed = 1)
#' s <- pe_split(x, test_frac = 0.5, seed = 1)
#' preds <- list(baseline_mean(x, s),
#'               baseline_linear(x, s, targets = x$sim$targets))
#' res <- pe_evaluate(preds, x, metrics = c("mse", "pds"))
#' pe_compare(res, reference = "mean", seed = 1)
pe_compare <- function(results, metrics = NULL, reference = NULL,
                       n_boot = 2000, conf = 0.95, p_adjust = "holm",
                       seed = NULL) {
  needed <- c("method", "perturbation", "metric", "value")
  if (!is.data.frame(results) || !all(needed %in% names(results))) {
    stop("`results` must be a data frame returned by pe_evaluate().",
         call. = FALSE)
  }
  if (anyDuplicated(results[c("method", "perturbation", "metric")]) > 0) {
    stop("`results` contains repeated method/perturbation/metric rows. ",
         "Compare methods evaluated once each on a single split.",
         call. = FALSE)
  }
  available <- unique(results$metric)
  if (is.null(metrics)) {
    metrics <- available
  }
  unknown <- setdiff(metrics, available)
  if (length(unknown) > 0) {
    stop("Metric(s) not found in `results`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  methods <- unique(results$method)
  if (length(methods) < 2) {
    stop("`results` must contain at least two methods.", call. = FALSE)
  }
  if (is.null(reference)) {
    pairs <- utils::combn(methods, 2, simplify = FALSE)
  } else {
    if (!is.character(reference) || length(reference) != 1 ||
        !reference %in% methods) {
      stop("`reference` must be one of the methods in `results`: ",
           paste(methods, collapse = ", "), call. = FALSE)
    }
    pairs <- lapply(setdiff(methods, reference), function(m) c(m, reference))
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) ||
      n_boot < 1) {
    stop("`n_boot` must be a single positive number.", call. = FALSE)
  }
  n_boot <- as.integer(n_boot)
  if (!is.numeric(conf) || length(conf) != 1 || is.na(conf) ||
      conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)

  rows <- with_seed(seed, do.call(rbind, lapply(metrics, function(met) {
    do.call(rbind, lapply(pairs, function(pr) {
      compare_pair(results, met, pr[1], pr[2], n_boot, conf)
    }))
  })))

  rows$p_adj <- stats::p.adjust(rows$p_value, method = p_adjust)

  higher_better <- metric_higher_is_better(rows$metric)
  significant <- !is.na(rows$ci_low) & (rows$ci_low > 0 | rows$ci_high < 0)
  a_better <- ifelse(higher_better, rows$diff > 0, rows$diff < 0)
  rows$favours <- ifelse(
    is.na(rows$diff) | rows$n < 2, NA_character_,
    ifelse(significant, ifelse(a_better, rows$method_a, rows$method_b),
           "neither")
  )
  rownames(rows) <- NULL
  rows
}

compare_pair <- function(results, metric, a, b, n_boot, conf) {
  sub <- results[results$metric == metric, , drop = FALSE]
  va <- stats::setNames(sub$value[sub$method == a],
                        sub$perturbation[sub$method == a])
  vb <- stats::setNames(sub$value[sub$method == b],
                        sub$perturbation[sub$method == b])
  common <- intersect(names(va), names(vb))
  d <- va[common] - vb[common]
  keep <- !is.na(d)
  common <- common[keep]
  d <- unname(d[keep])
  n <- length(d)
  ci <- boot_mean_ci(d, n_boot, conf)

  data.frame(
    metric = metric,
    method_a = a,
    method_b = b,
    n = n,
    mean_a = if (n > 0) mean(va[common]) else NA_real_,
    mean_b = if (n > 0) mean(vb[common]) else NA_real_,
    diff = if (n > 0) mean(d) else NA_real_,
    ci_low = ci[1],
    ci_high = ci[2],
    p_value = paired_p(d),
    stringsAsFactors = FALSE
  )
}

# Percentile bootstrap interval for the mean of `d`.
boot_mean_ci <- function(d, n_boot, conf) {
  n <- length(d)
  if (n < 2) {
    return(c(NA_real_, NA_real_))
  }
  idx <- matrix(sample.int(n, n * n_boot, replace = TRUE), nrow = n_boot)
  means <- rowMeans(matrix(d[idx], nrow = n_boot))
  alpha <- (1 - conf) / 2
  stats::quantile(means, c(alpha, 1 - alpha), names = FALSE)
}

# Two-sided paired Wilcoxon signed-rank test on the differences.
paired_p <- function(d) {
  if (length(d) < 2) {
    return(NA_real_)
  }
  if (all(d == 0)) {
    return(1)
  }
  suppressWarnings(stats::wilcox.test(d, mu = 0, exact = FALSE)$p.value)
}
