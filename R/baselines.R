#' Reference baselines for perturbation prediction
#'
#' Simple baselines that any useful model should beat. Reporting them
#' alongside a model is essential: recent benchmarks have shown that common
#' metrics can rank such baselines close to, or above, complex models.
#'
#' * `baseline_control()` predicts **no change**: every held-out perturbation
#'   is predicted to look like the average control cell.
#' * `baseline_mean()` predicts the **average perturbation response** seen in
#'   training: every held-out perturbation is predicted to look like the mean
#'   of the training perturbations' pseudobulk profiles (each perturbation
#'   weighted equally, regardless of its number of cells). It captures
#'   systematic differences shared by perturbed cells, but nothing
#'   perturbation-specific.
#' * `baseline_additive()` predicts a **combination** as the control profile
#'   plus the sum of its components' changes from control, each measured as a
#'   single perturbation in training. It needs an `"unseen_combo"` split.
#' * `baseline_linear()` is the **linear model** of Ahlmann-Eltze, Huber and
#'   Anders (Nat Methods 2025) for unseen single perturbations. Gene
#'   embeddings `G` are the top principal components of the centred training
#'   pseudobulk matrix `Y`; each perturbation is represented by the embedding
#'   of its target gene; and `Y` is modelled as
#'   `b + G W t(P)` with the ridge solution
#'   `W = (t(G) G + ridge I)^-1 t(G) (Y - b) P (t(P) P + ridge I)^-1`,
#'   where `b` is the row mean of `Y`. As `ridge` grows the prediction tends
#'   to `b`, the mean of the training single perturbations.
#'
#' `baseline_control()` and `baseline_mean()` make the same prediction for
#' every held-out perturbation, so they score at chance level (0.5 on
#' average) on the perturbation discrimination score of [pe_evaluate()].
#'
#' @param x A `pe_data` object.
#' @param split A `pe_split` object created from `x`.
#' @param n_comp For `baseline_linear()`: number of principal components used
#'   as gene embeddings (capped at the number of training perturbations
#'   minus one).
#' @param ridge For `baseline_linear()`: ridge penalty, a non-negative number.
#' @param targets For `baseline_linear()`: a named character vector mapping
#'   each perturbation label to its target gene (a row name of `x$expr`).
#'   By default the label itself is taken to be the target gene, as is usual
#'   for CRISPR screens labelled by gene symbol. For [pe_simulate()] data use
#'   `x$sim$targets`.
#'
#' @return A `pe_prediction` object (see [pe_prediction()]) whose `pred`
#'   element is a genes x held-out-perturbations matrix laid out like
#'   [pe_truth()].
#' @name baselines
#' @examples
#' x <- pe_simulate(seed = 1)
#'
#' # Unseen single perturbations
#' s <- pe_split(x, seed = 1)
#' preds <- list(
#'   baseline_control(x, s),
#'   baseline_mean(x, s),
#'   baseline_linear(x, s, targets = x$sim$targets)
#' )
#' pe_summary(pe_evaluate(preds, x, metrics = c("mse", "pearson_delta_pert")))
#'
#' # Unseen combinations
#' s2 <- pe_split(x, type = "unseen_combo", seed = 1)
#' pe_summary(pe_evaluate(baseline_additive(x, s2), x, metrics = "mse"))
NULL

#' @rdname baselines
#' @export
baseline_control <- function(x, split) {
  check_pe_data(x)
  check_pe_split(split, x)
  profile <- pe_pseudobulk(x, x$control)[, 1]
  pe_prediction(repeat_profile(profile, split$test), split, method = "control")
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
  pe_prediction(repeat_profile(profile, split$test), split, method = "mean")
}

#' @rdname baselines
#' @export
baseline_additive <- function(x, split) {
  check_pe_data(x)
  check_pe_split(split, x)
  sep <- x$combo_sep
  if (!all(is_combo(split$test, sep))) {
    stop("baseline_additive() predicts combinatorial perturbations only; ",
         "use a split made with type = \"unseen_combo\".", call. = FALSE)
  }
  parts <- combo_parts(split$test, sep)
  needed <- unique(unlist(parts))
  missing <- setdiff(needed, split$train)
  if (length(missing) > 0) {
    stop("These single perturbations are needed but not in the training ",
         "set: ", paste(utils::head(missing, 5), collapse = ", "),
         call. = FALSE)
  }
  ctrl <- pe_pseudobulk(x, x$control)[, 1]
  deltas <- pe_pseudobulk(x, needed) - ctrl
  pred <- vapply(parts, function(g) ctrl + rowSums(deltas[, g, drop = FALSE]),
                 numeric(length(ctrl)))
  pred <- matrix(pred, nrow = length(ctrl),
                 dimnames = list(names(ctrl), split$test))
  pe_prediction(pred, split, method = "additive")
}

#' @rdname baselines
#' @export
baseline_linear <- function(x, split, n_comp = 10, ridge = 0.1,
                            targets = NULL) {
  check_pe_data(x)
  check_pe_split(split, x)
  if (!is.numeric(n_comp) || length(n_comp) != 1 || is.na(n_comp) ||
      n_comp < 1) {
    stop("`n_comp` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(ridge) || length(ridge) != 1 || is.na(ridge) || ridge < 0) {
    stop("`ridge` must be a single non-negative number.", call. = FALSE)
  }
  sep <- x$combo_sep
  if (any(is_combo(split$test, sep))) {
    stop("baseline_linear() predicts single perturbations only; ",
         "use a split made with type = \"unseen_single\".", call. = FALSE)
  }
  train <- split$train[!is_combo(split$train, sep)]
  if (length(train) < 2) {
    stop("baseline_linear() needs at least two training single ",
         "perturbations.", call. = FALSE)
  }
  tg <- resolve_targets(c(train, split$test), targets, rownames(x$expr))

  y <- pe_pseudobulk(x, train)
  b <- rowMeans(y)
  y_centred <- y - b
  k <- min(as.integer(n_comp), ncol(y) - 1L)
  g <- svd(y_centred, nu = k, nv = 0)$u
  rownames(g) <- rownames(y)

  p_train <- g[tg[train], , drop = FALSE]
  p_test <- g[tg[split$test], , drop = FALSE]
  ident <- diag(k)
  w <- solve(crossprod(g) + ridge * ident, crossprod(g, y_centred)) %*%
    p_train %*% solve(crossprod(p_train) + ridge * ident)

  pred <- b + g %*% w %*% t(p_test)
  dimnames(pred) <- list(rownames(y), split$test)
  pe_prediction(pred, split, method = "linear")
}

# Helpers -------------------------------------------------------------------

repeat_profile <- function(profile, perturbations) {
  matrix(profile, nrow = length(profile), ncol = length(perturbations),
         dimnames = list(names(profile), perturbations))
}

resolve_targets <- function(perts, targets, genes) {
  if (is.null(targets)) {
    tg <- stats::setNames(perts, perts)
  } else {
    if (!is.character(targets) || is.null(names(targets))) {
      stop("`targets` must be a named character vector ",
           "(names = perturbation labels, values = target genes).",
           call. = FALSE)
    }
    missing <- setdiff(perts, names(targets))
    if (length(missing) > 0) {
      stop("`targets` has no entry for: ",
           paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
    }
    tg <- stats::setNames(unname(targets[perts]), perts)
  }
  unknown <- setdiff(tg, genes)
  if (length(unknown) > 0) {
    stop("Target genes not found among the measured genes: ",
         paste(utils::head(unknown, 5), collapse = ", "),
         ". Supply `targets` to map perturbation labels to genes.",
         call. = FALSE)
  }
  tg
}
