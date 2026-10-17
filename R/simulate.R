#' Simulate a small perturbation dataset
#'
#' Generates log-scale expression for control cells, single-gene
#' perturbations and additive combinations. Useful for examples, tests and
#' for checking how metrics behave when the ground truth is known.
#'
#' Each cell's expression is
#' `baseline + shared + effect[perturbation] + noise`, where `shared` is a
#' systematic shift applied to every perturbed (non-control) cell, and the
#' effect of a combination is the sum of its components' effects. Raising
#' `shared_effect` relative to `effect_sd` mimics datasets with strong
#' systematic variation, where the mean baseline looks deceptively good.
#'
#' Each single perturbation knocks down its own target gene (an effect of
#' `-2 * effect_sd`) in addition to changing `n_affected` random genes.
#' With `effect_spread > 0`, the whole effect of each single perturbation is
#' multiplied by its own strength, drawn from a log-normal distribution, so
#' that some perturbations are strong and others barely measurable, as in
#' real screens.
#'
#' @param n_genes Number of genes.
#' @param n_single Number of single-gene perturbations.
#' @param n_combo Number of two-gene combinations (at most
#'   `choose(n_single, 2)`).
#' @param cells_per_pert Cells per perturbation.
#' @param n_control Number of control cells.
#' @param n_affected Genes affected by each single perturbation.
#' @param effect_sd Standard deviation of perturbation-specific effects.
#' @param shared_effect Standard deviation of the systematic shift shared by
#'   all perturbed cells.
#' @param noise_sd Standard deviation of cell-level noise.
#' @param effect_spread Standard deviation of the log strength of single
#'   perturbations. The default 0 gives all perturbations the same strength.
#' @param seed Optional integer seed. The global RNG state is restored
#'   afterwards.
#'
#' @return A `pe_data` object with control label `"control"`, single
#'   perturbations named `P1`, `P2`, ... and combinations such as `"P1+P2"`.
#'   The simulated truth is stored in the extra element `sim`, a list with
#'   `baseline`, `shared`, the genes x perturbations matrix `effects`, and
#'   `targets`, a named vector giving each single perturbation's target gene
#'   (suitable for the `targets` argument of [baseline_linear()]), and
#'   `strength`, the strength multiplier of each single perturbation.
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' x
pe_simulate <- function(n_genes = 200, n_single = 20, n_combo = 10,
                        cells_per_pert = 30, n_control = 200,
                        n_affected = 20, effect_sd = 1,
                        shared_effect = 0.5, noise_sd = 0.5,
                        effect_spread = 0, seed = NULL) {
  if (n_single < 2) {
    stop("`n_single` must be at least 2.", call. = FALSE)
  }
  if (n_combo > choose(n_single, 2)) {
    stop("`n_combo` cannot exceed choose(n_single, 2).", call. = FALSE)
  }
  if (n_affected > n_genes) {
    stop("`n_affected` cannot exceed `n_genes`.", call. = FALSE)
  }
  if (n_single > n_genes) {
    stop("`n_single` cannot exceed `n_genes`.", call. = FALSE)
  }
  if (!is.numeric(effect_spread) || length(effect_spread) != 1 ||
      is.na(effect_spread) || effect_spread < 0) {
    stop("`effect_spread` must be a single non-negative number.",
         call. = FALSE)
  }
  with_seed(seed, simulate_impl(
    n_genes, n_single, n_combo, cells_per_pert, n_control,
    n_affected, effect_sd, shared_effect, noise_sd, effect_spread
  ))
}

simulate_impl <- function(n_genes, n_single, n_combo, cells_per_pert,
                          n_control, n_affected, effect_sd, shared_effect,
                          noise_sd, effect_spread = 0) {
  genes <- paste0("g", seq_len(n_genes))
  singles <- paste0("P", seq_len(n_single))

  baseline <- stats::rnorm(n_genes, mean = 2, sd = 1)
  shared <- stats::rnorm(n_genes, mean = 0, sd = shared_effect)

  # Each single perturbation knocks down one distinct target gene and changes
  # `n_affected` other randomly chosen genes.
  target_idx <- sample.int(n_genes, n_single)
  single_eff <- matrix(0, n_genes, n_single, dimnames = list(genes, singles))
  for (j in seq_len(n_single)) {
    hit <- sample.int(n_genes, n_affected)
    single_eff[hit, j] <- stats::rnorm(n_affected, mean = 0, sd = effect_sd)
    single_eff[target_idx[j], j] <- -2 * effect_sd
  }
  # Random draws for the strengths happen only when requested, so results
  # for effect_spread = 0 are identical to earlier versions.
  strength <- stats::setNames(rep(1, n_single), singles)
  if (effect_spread > 0) {
    strength[] <- exp(stats::rnorm(n_single, mean = 0, sd = effect_spread))
    single_eff <- sweep(single_eff, 2, strength, `*`)
  }

  combo_names <- character(0)
  combo_eff <- matrix(0, n_genes, 0, dimnames = list(genes, NULL))
  if (n_combo > 0) {
    pairs <- utils::combn(n_single, 2)
    chosen <- pairs[, sample.int(ncol(pairs), n_combo), drop = FALSE]
    combo_names <- paste(singles[chosen[1, ]], singles[chosen[2, ]], sep = "+")
    combo_eff <- single_eff[, chosen[1, ], drop = FALSE] +
      single_eff[, chosen[2, ], drop = FALSE]
    colnames(combo_eff) <- combo_names
  }
  effects <- cbind(single_eff, combo_eff)

  labels <- c(rep("control", n_control),
              rep(colnames(effects), each = cells_per_pert))
  mean_profile <- cbind(
    matrix(baseline, n_genes, n_control),
    (baseline + shared + effects)[, rep(seq_len(ncol(effects)),
                                        each = cells_per_pert), drop = FALSE]
  )
  noise <- matrix(stats::rnorm(length(mean_profile), sd = noise_sd),
                  nrow = n_genes)
  expr <- mean_profile + noise
  dimnames(expr) <- list(genes, paste0("cell", seq_along(labels)))

  out <- pe_data(expr, labels, control = "control", combo_sep = "+")
  out$sim <- list(baseline = stats::setNames(baseline, genes),
                  shared = stats::setNames(shared, genes),
                  effects = effects,
                  targets = stats::setNames(genes[target_idx], singles),
                  strength = strength)
  out
}
