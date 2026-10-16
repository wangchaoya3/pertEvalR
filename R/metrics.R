#' Evaluate perturbation predictions
#'
#' Scores one or more predictions against the observed pseudobulk profiles of
#' the held-out perturbations, returning one value per perturbation and metric.
#'
#' All metrics compare pseudobulk (mean) profiles. Let `pred` and `truth` be
#' the predicted and observed profiles of a held-out perturbation, `ctrl` the
#' mean control profile and `centroid` the mean of the training
#' perturbations' profiles (the prediction of [baseline_mean()]).
#'
#' * `"mse"`: mean squared error over all genes.
#' * `"mse_top"`: mean squared error over the `top_n` genes that change most
#'   in the observed data (largest `|truth - ctrl|`), chosen per perturbation.
#' * `"pearson_delta"`: Pearson correlation of `pred - ctrl` with
#'   `truth - ctrl`, over all genes. Rewards predicting the direction of
#'   change from control, including changes shared by all perturbations.
#' * `"pearson_delta_top"`: as above, over the `top_n` genes.
#' * `"pearson_delta_pert"`: Pearson correlation of `pred - centroid` with
#'   `truth - centroid`. Using the perturbed centroid as reference removes
#'   the systematic variation shared by perturbed cells, so only
#'   perturbation-specific effects count (the idea behind Systema, Vinas Torne
#'   et al., Nat Biotechnol 2025).
#' * `"pds"`: perturbation discrimination score. For each held-out
#'   perturbation, the fraction of the other held-out perturbations whose
#'   observed change from control is further from the predicted change than
#'   its own observed change is (ties count one half). 1 is perfect, 0.5 is
#'   chance. Needs at least two held-out perturbations.
#'
#' Correlation metrics are `NA` when either vector is constant. In particular
#' `baseline_control()` has no `pearson_delta` and `baseline_mean()` has no
#' `pearson_delta_pert`: those baselines predict no change from the reference
#' by construction.
#'
#' @param prediction A `pe_prediction` object, or a list of them.
#' @param x The `pe_data` object the predictions were made for.
#' @param metrics Metrics to compute; see Details. Defaults to all.
#' @param top_n Number of most-changed genes used by the `_top` metrics.
#' @param distance Distance used by `"pds"`: `"l1"`, `"l2"` or `"cosine"`
#'   (one minus cosine similarity).
#'
#' @return A data frame with columns `method`, `perturbation`, `metric` and
#'   `value`, one row per prediction, held-out perturbation and metric. Use
#'   [pe_summary()] to average over perturbations.
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' res <- pe_evaluate(list(baseline_control(x, s), baseline_mean(x, s)), x)
#' head(res)
#' pe_summary(res)
pe_evaluate <- function(prediction, x,
                        metrics = c("mse", "mse_top", "pearson_delta",
                                    "pearson_delta_top", "pearson_delta_pert",
                                    "pds"),
                        top_n = 20,
                        distance = c("l1", "l2", "cosine")) {
  check_pe_data(x)
  if (inherits(prediction, "pe_prediction")) {
    prediction <- list(prediction)
  }
  is_pred <- is.list(prediction) && length(prediction) > 0 &&
    all(vapply(prediction, inherits, logical(1), what = "pe_prediction"))
  if (!is_pred) {
    stop("`prediction` must be a `pe_prediction` object or a list of them.",
         call. = FALSE)
  }
  metrics <- match.arg(metrics, several.ok = TRUE)
  distance <- match.arg(distance)
  if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) || top_n < 1) {
    stop("`top_n` must be a single positive number.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  results <- lapply(prediction, evaluate_one, x = x, metrics = metrics,
                    top_n = top_n, distance = distance)
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  attr(out, "settings") <- list(top_n = top_n, distance = distance)
  out
}

evaluate_one <- function(p, x, metrics, top_n, distance) {
  split <- p$split
  check_pe_split(split, x)
  truth <- pe_truth(x, split)
  missing <- setdiff(rownames(truth), rownames(p$pred))
  if (length(missing) > 0) {
    stop("Prediction \"", p$method, "\" lacks ", length(missing),
         " gene(s) present in the data, e.g. ",
         paste(utils::head(missing, 3), collapse = ", "), ".", call. = FALSE)
  }
  pred <- p$pred[rownames(truth), colnames(truth), drop = FALSE]
  ctrl <- pe_pseudobulk(x, x$control)[, 1]
  centroid <- NULL
  if (length(split$train) > 0) {
    centroid <- rowMeans(pe_pseudobulk(x, split$train))
  }
  top <- top_genes(truth, ctrl, min(top_n, nrow(truth)))

  values <- lapply(metrics, function(m) {
    switch(m,
      mse = colMeans((pred - truth)^2),
      mse_top = vapply(seq_len(ncol(truth)), function(j) {
        g <- top[[j]]
        mean((pred[g, j] - truth[g, j])^2)
      }, numeric(1)),
      pearson_delta = col_cor(pred - ctrl, truth - ctrl),
      pearson_delta_top = col_cor(pred - ctrl, truth - ctrl, top),
      pearson_delta_pert = if (is.null(centroid)) {
        rep(NA_real_, ncol(truth))
      } else {
        col_cor(pred - centroid, truth - centroid)
      },
      pds = pds_scores(pred - ctrl, truth - ctrl, distance)
    )
  })

  data.frame(
    method = p$method,
    perturbation = rep(colnames(truth), times = length(metrics)),
    metric = rep(metrics, each = ncol(truth)),
    value = unlist(lapply(values, unname), use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

# Indices of the `n` genes with the largest observed change, per column.
top_genes <- function(truth, ctrl, n) {
  lapply(seq_len(ncol(truth)), function(j) {
    order(abs(truth[, j] - ctrl), decreasing = TRUE)[seq_len(n)]
  })
}

# Column-wise Pearson correlation; optional per-column row subsets.
col_cor <- function(a, b, rows = NULL) {
  vapply(seq_len(ncol(a)), function(j) {
    g <- if (is.null(rows)) seq_len(nrow(a)) else rows[[j]]
    safe_cor(a[g, j], b[g, j])
  }, numeric(1))
}

safe_cor <- function(a, b) {
  if (length(a) < 2 || stats::sd(a) == 0 || stats::sd(b) == 0) {
    return(NA_real_)
  }
  stats::cor(a, b)
}

# Perturbation discrimination score from predicted and observed deltas.
pds_scores <- function(pred_delta, truth_delta, distance) {
  n <- ncol(truth_delta)
  if (n < 2) {
    return(rep(NA_real_, n))
  }
  dist <- vapply(seq_len(n), function(i) {
    delta_distance(pred_delta[, i], truth_delta, distance)
  }, numeric(n))
  # dist[j, i] = distance between prediction i and observation j
  vapply(seq_len(n), function(i) {
    own <- dist[i, i]
    others <- dist[-i, i]
    mean(others > own) + 0.5 * mean(others == own)
  }, numeric(1))
}

# Distance from one vector to every column of a matrix.
delta_distance <- function(v, m, distance) {
  switch(distance,
    l1 = colSums(abs(m - v)),
    l2 = sqrt(colSums((m - v)^2)),
    cosine = {
      nv <- sqrt(sum(v^2))
      nm <- sqrt(colSums(m^2))
      sim <- colSums(m * v) / (nv * nm)
      sim[!is.finite(sim)] <- 0
      1 - sim
    }
  )
}

#' Summarise evaluation results
#'
#' Averages the per-perturbation results of [pe_evaluate()] for each method
#' and metric.
#'
#' @param results A data frame returned by [pe_evaluate()].
#'
#' @return A data frame with columns `method`, `metric`, `mean`, `sd`, `n`
#'   (perturbations with a value) and `n_na` (perturbations where the metric
#'   is undefined).
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' pe_summary(pe_evaluate(baseline_mean(x, s), x))
pe_summary <- function(results) {
  needed <- c("method", "metric", "value")
  if (!is.data.frame(results) || !all(needed %in% names(results))) {
    stop("`results` must be a data frame returned by pe_evaluate().",
         call. = FALSE)
  }
  keys <- unique(results[c("method", "metric")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    v <- results$value[results$method == keys$method[i] &
                         results$metric == keys$metric[i]]
    ok <- !is.na(v)
    data.frame(
      method = keys$method[i],
      metric = keys$metric[i],
      mean = if (any(ok)) mean(v[ok]) else NA_real_,
      sd = if (sum(ok) > 1) stats::sd(v[ok]) else NA_real_,
      n = sum(ok),
      n_na = sum(!ok),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
