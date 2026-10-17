#' How well can each metric separate informative from uninformative predictions?
#'
#' Summarises a [pe_bounds()] table into two meta-metrics per metric, which
#' describe the metric on a given dataset rather than any model:
#'
#' * **Dynamic range fraction** (DRF), following Miller et al. (2026): the
#'   share of the metric's ideal range that lies between the lower and the
#'   upper bound,
#'   `drf = (mean(upper) - mean(lower)) / (ideal - mean(lower))`,
#'   where `ideal` is the best possible value (0 for `mse` and `mse_top`, 1
#'   otherwise). Values near 0 mean that even a technical replicate barely
#'   improves on an uninformative prediction, so the metric cannot reward a
#'   good model on these data.
#' * **Bound discrimination score** (BDS), following Vollenweider and
#'   Buehlmann (2026): the fraction of perturbations for which the metric
#'   scores the technical replicate better than the best uninformative
#'   baseline (ties count one half). 0.5 means the metric cannot tell the two
#'   apart; 1 means it always ranks the replicate first.
#'
#' Both are computed from the bounds of `pe_bounds()`: the lower bound is the
#' better of the no-change and training-mean baselines, the upper bound the
#' split-half replicate corrected to full sample size. Values therefore differ
#' in detail from those of the original publications, which construct their
#' controls differently.
#'
#' @param bounds A data frame returned by [pe_bounds()], optionally restricted
#'   with [pe_subset()].
#'
#' @return A data frame with one row per metric and columns `metric`, `n`
#'   (perturbations with both bounds), `lower`, `upper` (means), `ideal`,
#'   `drf` and `bds`.
#' @seealso [pe_bounds()], [pe_calibrate()]
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' pe_meta_metrics(pe_bounds(x, s, seed = 1))
pe_meta_metrics <- function(bounds) {
  needed <- c("perturbation", "metric", "lower", "upper")
  if (!is.data.frame(bounds) || !all(needed %in% names(bounds))) {
    stop("`bounds` must be a data frame returned by pe_bounds().",
         call. = FALSE)
  }
  metrics <- unique(bounds$metric)
  rows <- lapply(metrics, function(m) {
    d <- bounds[bounds$metric == m & !is.na(bounds$lower) &
                  !is.na(bounds$upper), , drop = FALSE]
    hb <- metric_higher_is_better(m)
    ideal <- if (hb) 1 else 0
    lower <- upper <- drf <- bds <- NA_real_
    if (nrow(d) > 0) {
      lower <- mean(d$lower)
      upper <- mean(d$upper)
      if (ideal != lower) {
        drf <- (upper - lower) / (ideal - lower)
      }
      better <- if (hb) d$upper > d$lower else d$upper < d$lower
      bds <- mean(better + 0.5 * (d$upper == d$lower))
    }
    data.frame(metric = m, n = nrow(d), lower = lower, upper = upper,
               ideal = ideal, drf = drf, bds = bds,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
