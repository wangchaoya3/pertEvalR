#' Reference baselines for perturbation prediction
#'
#' Uninformative baselines that any useful model should beat. Reporting them
#' alongside a model is essential: recent benchmarks have shown that common
#' metrics can rank these baselines close to, or above, complex models.
#'
#' * `baseline_control()` predicts **no change**: every held-out perturbation
#'   is predicted to look like the average control cell.
#' * `baseline_mean()` predicts the **average perturbation response** seen in
#'   training: every held-out perturbation is predicted to look like the mean
#'   of the training perturbations' pseudobulk profiles (each perturbation
#'   weighted equally, regardless of its number of cells). It captures
#'   systematic differences shared by perturbed cells, but nothing
#'   perturbation-specific.
#'
#' Neither baseline varies between held-out perturbations, so both score
#' zero on any metric that asks whether predictions tell perturbations apart.
#'
#' @param x A `pe_data` object.
#' @param split A `pe_split` object created from `x`.
#'
#' @return A `pe_prediction` object: a list with `pred` (a genes x
#'   test-perturbations matrix, with the same layout as [pe_truth()]),
#'   `method` and `split`.
#' @name baselines
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' p_ctrl <- baseline_control(x, s)
#' p_mean <- baseline_mean(x, s)
#' truth <- pe_truth(x, s)
#'
#' # Mean squared error of each baseline on the held-out perturbations
#' mean((p_ctrl$pred - truth)^2)
#' mean((p_mean$pred - truth)^2)
NULL

#' @rdname baselines
#' @export
baseline_control <- function(x, split) {
  check_pe_data(x)
  check_pe_split(split, x)
  profile <- pe_pseudobulk(x, x$control)[, 1]
  new_pe_prediction(repeat_profile(profile, split$test),
                    method = "control", split = split)
}

#' @rdname baselines
#' @export
baseline_mean <- function(x, split) {
  check_pe_data(x)
  check_pe_split(split, x)
  if (length(split$train) == 0) {
    stop("The split has no training perturbations.", call. = FALSE)
  }
  profile <- rowMeans(pe_pseudobulk(x, split$train))
  new_pe_prediction(repeat_profile(profile, split$test),
                    method = "mean", split = split)
}

# Prediction objects --------------------------------------------------------

repeat_profile <- function(profile, perturbations) {
  matrix(profile, nrow = length(profile), ncol = length(perturbations),
         dimnames = list(names(profile), perturbations))
}

new_pe_prediction <- function(pred, method, split) {
  stopifnot(is.matrix(pred), identical(colnames(pred), split$test))
  structure(list(pred = pred, method = method, split = split),
            class = "pe_prediction")
}

#' @export
print.pe_prediction <- function(x, ...) {
  cat("<pe_prediction: ", x$method, ">\n", sep = "")
  cat("  ", nrow(x$pred), " genes x ", ncol(x$pred),
      " held-out perturbations (", x$split$type, " split)\n", sep = "")
  invisible(x)
}
