#' Diagnose the perturbation signal in a dataset
#'
#' Describes how much perturbation-specific signal a dataset carries before
#' any model is evaluated on it. The key quantity is systematic variation:
#' when most perturbations shift expression in the same direction, metrics
#' that compare predictions with the change from control reward predicting
#' that common shift (see [baseline_mean()]), whatever the perturbation.
#'
#' For each perturbation, let `delta` be its pseudobulk profile minus the
#' control profile, and `average` the mean of `delta` over all perturbations
#' in the dataset.
#'
#' * `cosine_to_average`: cosine similarity between `delta` and `average`.
#'   Values near 1 mean the perturbation mostly reproduces the common shift.
#' * `effect_size`: root mean square of `delta` over genes (the typical
#'   change per gene).
#' * `effective_genes`: participation ratio
#'   sum(delta^2)^2 / sum(delta^4), roughly the number of genes that carry
#'   the change. Small values mean the effect is concentrated in a few genes.
#'   Sampling noise spreads over all genes and inflates this number when the
#'   effect is weak, so compare it between perturbations or datasets rather
#'   than reading it as a count of affected genes.
#' * `reliability`: the cells are split at random into two halves and the
#'   two half-sample deltas are correlated over genes; the split-half
#'   correlation r is converted to full-sample reliability with the
#'   Spearman-Brown formula 2r / (1 + r) (left unchanged when r is zero or
#'   negative). Low values mean the change is hard to tell apart from
#'   sampling noise. `NA` for perturbations with fewer than two cells.
#'   Reliability grows with the number of cells, so it describes how well an
#'   effect is measured in this dataset, not how strong it is biologically.
#'   [pe_subset()] restricts an evaluation to reliable perturbations.
#'
#' The dataset summary reports the mean and median cosine and the
#' systematic fraction: the share of the total squared change of the
#' selected perturbations that lies along the direction of `average`.
#'
#' @param x A `pe_data` object with at least two perturbations.
#' @param perturbations Perturbations to describe; defaults to all. The
#'   average shift is always computed over all perturbations in `x`.
#' @param seed Optional integer seed for the random halving used by
#'   `reliability`. The global random number generator state is restored
#'   afterwards.
#'
#' @return An object of class `pe_signal`: a list with `per_perturbation`, a
#'   data frame with columns `perturbation`, `n_cells`, `effect_size`,
#'   `cosine_to_average`, `effective_genes` and `reliability`, and
#'   `summary`, a list with `n_perturbations`, `mean_cosine`,
#'   `median_cosine`, `systematic_fraction` and `median_reliability`.
#' @export
#' @examples
#' strong <- pe_simulate(shared_effect = 2, seed = 1)
#' weak <- pe_simulate(shared_effect = 0, seed = 1)
#' pe_signal(strong, seed = 1)
#' pe_signal(weak, seed = 1)
pe_signal <- function(x, perturbations = NULL, seed = NULL) {
  check_pe_data(x)
  all_p <- pe_perturbations(x)
  if (length(all_p) < 2) {
    stop("`x` must contain at least two perturbations.", call. = FALSE)
  }
  if (is.null(perturbations)) {
    perturbations <- all_p
  }
  perturbations <- as.character(perturbations)
  unknown <- setdiff(perturbations, all_p)
  if (length(unknown) > 0) {
    stop("Unknown perturbation(s): ",
         paste(utils::head(unknown, 5), collapse = ", "), call. = FALSE)
  }

  pb <- pe_pseudobulk(x, c(x$control, all_p))
  ctrl <- pb[, 1]
  delta_all <- pb[, -1, drop = FALSE] - ctrl
  avg <- rowMeans(delta_all)
  d <- delta_all[, perturbations, drop = FALSE]

  ss <- colSums(d^2)
  norm_d <- sqrt(ss)
  norm_avg <- sqrt(sum(avg^2))
  cosine <- colSums(d * avg) / (norm_d * norm_avg)
  cosine[!is.finite(cosine)] <- NA_real_
  eff_genes <- ss^2 / colSums(d^4)
  eff_genes[!is.finite(eff_genes)] <- NA_real_

  reliability <- with_seed(seed, split_half_reliability(x, perturbations,
                                                        ctrl))
  n_cells <- as.vector(table(factor(x$perturbation, levels = perturbations)))

  per <- data.frame(
    perturbation = perturbations,
    n_cells = n_cells,
    effect_size = unname(norm_d / sqrt(nrow(d))),
    cosine_to_average = unname(cosine),
    effective_genes = unname(eff_genes),
    reliability = unname(reliability),
    stringsAsFactors = FALSE
  )

  proj <- colSums(d * avg) / norm_avg
  summary <- list(
    n_perturbations = length(perturbations),
    mean_cosine = mean(per$cosine_to_average, na.rm = TRUE),
    median_cosine = stats::median(per$cosine_to_average, na.rm = TRUE),
    systematic_fraction = if (norm_avg > 0 && sum(ss) > 0) {
      sum(proj^2) / sum(ss)
    } else {
      NA_real_
    },
    median_reliability = if (any(!is.na(per$reliability))) {
      stats::median(per$reliability, na.rm = TRUE)
    } else {
      NA_real_
    }
  )
  structure(list(per_perturbation = per, summary = summary),
            class = "pe_signal")
}

split_half_reliability <- function(x, perts, ctrl) {
  suffix <- "__pe_signal_half"
  x_half <- split_cells_in_half(x, perts, suffix)
  out <- stats::setNames(rep(NA_real_, length(perts)), perts)
  halved <- perts[paste0(perts, suffix) %in% x_half$perturbation]
  if (length(halved) == 0) {
    return(out)
  }
  a <- pe_pseudobulk(x_half, paste0(halved, suffix)) - ctrl
  b <- pe_pseudobulk(x_half, halved) - ctrl
  r <- vapply(seq_along(halved), function(j) safe_cor(a[, j], b[, j]),
              numeric(1))
  rel <- ifelse(!is.na(r) & r > 0, 2 * r / (1 + r), r)
  out[halved] <- rel
  out
}

#' @export
print.pe_signal <- function(x, ...) {
  s <- x$summary
  cat("<pe_signal>\n")
  cat("  perturbations:        ", s$n_perturbations, "\n", sep = "")
  cat(sprintf("  cosine to average:     mean %.2f, median %.2f\n",
              s$mean_cosine, s$median_cosine))
  cat(sprintf("  systematic fraction:   %.2f of the squared change lies along the average shift\n",
              s$systematic_fraction))
  cat(sprintf("  median reliability:    %.2f\n", s$median_reliability))
  invisible(x)
}
