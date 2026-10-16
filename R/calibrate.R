#' Lower and upper reference bounds for each metric
#'
#' Every metric can be flattered by some uninformative prediction: correlation
#' with the change from control rewards the training mean, and correlation
#' around the perturbed centroid can reward predicting no change. `pe_bounds()`
#' measures, for each held-out perturbation and metric, how well such
#' predictions already score (the lower bound) and how well an independent
#' measurement of the same perturbation scores (the upper bound), so that
#' model scores can be placed between the two with [pe_calibrate()].
#'
#' * **Lower bound**: for each metric, whichever of [baseline_control()] and
#'   [baseline_mean()] scores better on average over the held-out
#'   perturbations (higher mean, or lower for `mse` and `mse_top`); its
#'   per-perturbation values are the lower bound and `lower_method` records
#'   which baseline it is. Baselines whose value is undefined are ignored. The
#'   choice is made once per metric rather than per perturbation: taking the
#'   better of two noisy scores separately for every perturbation would bias
#'   the bound upwards.
#' * **Upper bound**: a technical replicate. The cells of each held-out
#'   perturbation are split at random into two halves and the pseudobulk of one
#'   half is scored as a prediction of the other half (`upper_raw`). Each half
#'   holds only half of the cells, so this raw score is noisier than a
#'   comparison with all cells; `upper` corrects for that:
#'   * `mse`, `mse_top`: the squared difference of two half-sample means has
#'     four times the sampling variance of one full-sample mean, so
#'     `upper = upper_raw / 4` (for `mse_top` approximately, as its genes are
#'     chosen from the data).
#'   * `pearson_delta`, `pearson_delta_top`, `pearson_delta_pert`: with
#'     split-half correlation r, the Spearman-Brown reliability of the
#'     full-sample profile is 2r / (1 + r), and the best correlation any
#'     prediction can reach with it is the square root of that,
#'     upper = sqrt(2r / (1 + r)) (left unchanged when r is zero or
#'     negative).
#'   * `pds`: no simple correction exists; `upper = upper_raw`.
#'
#'   The corrections assume independent cells and equal halves; scores
#'   slightly above the upper bound are therefore possible.
#'
#' This follows the bounds-based calibration proposed in the "Signal, Bounds
#' and Baselines" principles for evaluating perturbation models (2026).
#'
#' @inheritParams pe_evaluate
#' @param x A `pe_data` object.
#' @param split A `pe_split` object created from `x`. Every held-out
#'   perturbation needs at least two cells.
#' @param seed Optional integer seed for the random halving of cells. The
#'   global random number generator state is restored afterwards.
#'
#' @return A data frame with columns `perturbation`, `metric`, `lower`,
#'   `lower_method`, `upper` and `upper_raw`, one row per held-out
#'   perturbation and metric. Its `"settings"` attribute records `top_n` and
#'   `distance`.
#' @seealso [pe_calibrate()]
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' b <- pe_bounds(x, s, metrics = c("mse", "pearson_delta"), seed = 1)
#' head(b)
pe_bounds <- function(x, split,
                      metrics = c("mse", "mse_top", "pearson_delta",
                                  "pearson_delta_top", "pearson_delta_pert",
                                  "pds"),
                      top_n = 20, distance = c("l1", "l2", "cosine"),
                      seed = NULL) {
  check_pe_data(x)
  check_pe_split(split, x)
  metrics <- match.arg(metrics, several.ok = TRUE)
  distance <- match.arg(distance)

  n_cells <- table(factor(x$perturbation, levels = split$test))
  too_few <- names(n_cells)[n_cells < 2]
  if (length(too_few) > 0) {
    stop("These held-out perturbations have fewer than two cells, so no ",
         "replicate bound can be formed: ",
         paste(utils::head(too_few, 5), collapse = ", "), call. = FALSE)
  }

  # Lower bound: the uninformative baseline that is better on average, chosen
  # once per metric
  base <- pe_evaluate(list(baseline_control(x, split), baseline_mean(x, split)),
                      x, metrics = metrics, top_n = top_n, distance = distance)
  lower <- do.call(rbind, lapply(metrics, function(m) {
    d <- base[base$metric == m, ]
    avg <- vapply(split(d$value, d$method),
                  function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE),
                  numeric(1))
    avg <- avg[!is.na(avg)]
    if (length(avg) == 0) {
      return(data.frame(perturbation = split$test, metric = m,
                        lower = NA_real_, lower_method = NA_character_,
                        stringsAsFactors = FALSE))
    }
    best <- names(avg)[if (metric_higher_is_better(m)) which.max(avg)
                       else which.min(avg)]
    db <- d[d$method == best, ]
    data.frame(perturbation = db$perturbation, metric = m, lower = db$value,
               lower_method = best, stringsAsFactors = FALSE)
  }))

  # Upper bound: split-half technical replicate, corrected to full sample size
  upper <- with_seed(seed, replicate_scores(x, split, metrics, top_n,
                                            distance))
  upper$upper_raw <- upper$value
  upper$upper <- full_sample_bound(upper$value, upper$metric)

  out <- merge(lower, upper[c("perturbation", "metric", "upper", "upper_raw")],
               by = c("perturbation", "metric"), all = TRUE, sort = FALSE)
  out <- out[order(match(out$metric, metrics),
                   match(out$perturbation, split$test)), ]
  rownames(out) <- NULL
  attr(out, "settings") <- list(top_n = as.integer(top_n),
                                distance = distance)
  out
}

# Convert a split-half replicate score into a bound for comparisons with the
# full sample (see the pe_bounds() documentation).
full_sample_bound <- function(value, metric) {
  out <- value
  is_mse <- metric %in% c("mse", "mse_top")
  out[is_mse] <- value[is_mse] / 4
  is_cor <- metric %in% c("pearson_delta", "pearson_delta_top",
                          "pearson_delta_pert")
  pos <- is_cor & !is.na(value) & value > 0
  out[pos] <- sqrt(2 * value[pos] / (1 + value[pos]))
  out
}

# Score one random half of each held-out perturbation's cells as a
# prediction of the other half.
replicate_scores <- function(x, split, metrics, top_n, distance) {
  suffix <- "__pe_replicate_half"
  if (any(endsWith(unique(x$perturbation), suffix))) {
    stop("Perturbation labels must not end with \"", suffix, "\".",
         call. = FALSE)
  }
  labels <- x$perturbation
  for (p in split$test) {
    idx <- which(labels == p)
    half <- idx[sample.int(length(idx), floor(length(idx) / 2))]
    labels[half] <- paste0(p, suffix)
  }
  x_half <- x
  x_half$perturbation <- labels

  pred <- pe_pseudobulk(x_half, paste0(split$test, suffix))
  colnames(pred) <- split$test
  rep_pred <- pe_prediction(pred, split, method = "replicate")
  pe_evaluate(rep_pred, x_half, metrics = metrics, top_n = top_n,
              distance = distance)
}

#' Calibrate scores between reference bounds
#'
#' Places each method's score between the lower bound (best uninformative
#' baseline) and the upper bound (technical replicate) from [pe_bounds()]:
#'
#' `calibrated = (score - lower) / (upper - lower)`
#'
#' so that 0 means "no better than an uninformative baseline" and 1 means "as
#' good as an independent measurement of the same perturbation". Negative
#' values are worse than the best uninformative baseline. The same formula
#' applies to metrics where lower is better (`mse`, `mse_top`), because both
#' the numerator and the denominator change sign.
#'
#' `calibrated` uses the means over perturbations of the score and of the two
#' bounds, which is stable. `calibrated_median` is the median of
#' per-perturbation calibrated scores, which is more robust to a few
#' perturbations with extreme values but undefined where a perturbation's
#' bounds coincide. Calibration is `NA` when the upper bound is not better
#' than the lower bound, because the metric then has no room to show
#' improvement on this data.
#'
#' @param results A data frame returned by [pe_evaluate()].
#' @param bounds A data frame returned by [pe_bounds()] for the same data,
#'   split and evaluation settings.
#'
#' @return A data frame with one row per method and metric and columns
#'   `method`, `metric`, `n` (perturbations used), `score`, `lower`, `upper`
#'   (means over those perturbations), `calibrated` and `calibrated_median`.
#' @seealso [pe_bounds()]
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' res <- pe_evaluate(baseline_linear(x, s, targets = x$sim$targets), x)
#' b <- pe_bounds(x, s, seed = 1)
#' pe_calibrate(res, b)
pe_calibrate <- function(results, bounds) {
  needed <- c("method", "perturbation", "metric", "value")
  if (!is.data.frame(results) || !all(needed %in% names(results))) {
    stop("`results` must be a data frame returned by pe_evaluate().",
         call. = FALSE)
  }
  if (!is.data.frame(bounds) ||
      !all(c("perturbation", "metric", "lower", "upper") %in% names(bounds))) {
    stop("`bounds` must be a data frame returned by pe_bounds().",
         call. = FALSE)
  }
  s_res <- attr(results, "settings")
  s_bnd <- attr(bounds, "settings")
  if (!is.null(s_res) && !is.null(s_bnd) &&
      (!identical(as.integer(s_res$top_n), as.integer(s_bnd$top_n)) ||
       !identical(s_res$distance, s_bnd$distance))) {
    stop("`results` and `bounds` were computed with different `top_n` or ",
         "`distance` settings.", call. = FALSE)
  }
  if (!setequal(unique(results$perturbation), unique(bounds$perturbation))) {
    stop("`results` and `bounds` cover different perturbations; compute ",
         "them on the same split.", call. = FALSE)
  }
  missing_metrics <- setdiff(unique(results$metric), unique(bounds$metric))
  if (length(missing_metrics) > 0) {
    stop("`bounds` lacks metric(s): ",
         paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  m <- merge(results[needed], bounds[c("perturbation", "metric", "lower",
                                       "upper")],
             by = c("perturbation", "metric"), all.x = TRUE, sort = FALSE)

  keys <- unique(results[c("method", "metric")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    d <- m[m$method == keys$method[i] & m$metric == keys$metric[i], ]
    d <- d[!is.na(d$value) & !is.na(d$lower) & !is.na(d$upper), ]
    hb <- metric_higher_is_better(keys$metric[i])
    score <- lower <- upper <- cal <- cal_med <- NA_real_
    if (nrow(d) > 0) {
      score <- mean(d$value)
      lower <- mean(d$lower)
      upper <- mean(d$upper)
      cal <- calibrate_value(score, lower, upper, hb)
      per <- calibrate_value(d$value, d$lower, d$upper, hb)
      if (any(!is.na(per))) cal_med <- stats::median(per, na.rm = TRUE)
    }
    data.frame(method = keys$method[i], metric = keys$metric[i],
               n = nrow(d), score = score, lower = lower, upper = upper,
               calibrated = cal, calibrated_median = cal_med,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# (score - lower) / (upper - lower), NA where the upper bound is not better
# than the lower bound in the metric's direction.
calibrate_value <- function(score, lower, upper, higher_better) {
  headroom <- if (higher_better) upper > lower else upper < lower
  ifelse(headroom, (score - lower) / (upper - lower), NA_real_)
}
