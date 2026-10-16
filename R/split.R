#' Split perturbations into training and held-out sets
#'
#' Defines which perturbations a model may learn from and which it must
#' predict. Splits are made at the level of perturbations, never cells, so
#' that a held-out perturbation contributes no information to training.
#'
#' Two split types are available:
#'
#' * `"unseen_single"`: a fraction of the single-gene perturbations is held
#'   out. Any combinatorial perturbation that involves a held-out gene is
#'   removed from the training set as well, because it would otherwise leak
#'   information about that gene.
#' * `"unseen_combo"`: a fraction of the combinatorial perturbations is held
#'   out. Only combinations whose individual genes were all measured as single
#'   perturbations are eligible, so the task is to predict how known effects
#'   combine.
#'
#' The control label is never part of `train` or `test`; baselines and models
#' always have access to the control cells.
#'
#' @param x A `pe_data` object.
#' @param type The split type, `"unseen_single"` or `"unseen_combo"`.
#' @param test_frac Fraction of eligible perturbations to hold out. At least
#'   one perturbation is always held out.
#' @param seed Optional integer seed for a reproducible split. The global
#'   random number generator state is restored afterwards.
#'
#' @return An object of class `pe_split` with elements `train`, `test`,
#'   `control`, `type`, `test_frac` and `seed`.
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' pe_split(x, type = "unseen_single", seed = 42)
#' pe_split(x, type = "unseen_combo", seed = 42)
pe_split <- function(x, type = c("unseen_single", "unseen_combo"),
                     test_frac = 0.2, seed = NULL) {
  check_pe_data(x)
  type <- match.arg(type)
  if (!is.numeric(test_frac) || length(test_frac) != 1 ||
      is.na(test_frac) || test_frac <= 0 || test_frac >= 1) {
    stop("`test_frac` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  perts <- pe_perturbations(x, type = "all")
  singles <- pe_perturbations(x, type = "single")
  combos <- pe_perturbations(x, type = "combo")

  if (type == "unseen_single") {
    eligible <- singles
    if (length(eligible) < 2) {
      stop("At least two single perturbations are needed for an ",
           "\"unseen_single\" split.", call. = FALSE)
    }
  } else {
    if (length(combos) == 0) {
      stop("No combinatorial perturbations found (labels containing \"",
           x$combo_sep, "\").", call. = FALSE)
    }
    has_all_parts <- vapply(combo_parts(combos, x$combo_sep),
                            function(g) all(g %in% singles), logical(1))
    eligible <- combos[has_all_parts]
    if (length(eligible) == 0) {
      stop("No combination has all of its genes measured as single ",
           "perturbations, so an \"unseen_combo\" split is not possible.",
           call. = FALSE)
    }
  }

  n_test <- max(1L, as.integer(round(test_frac * length(eligible))))
  if (type == "unseen_single") {
    # Keep at least one single perturbation for training.
    n_test <- min(n_test, length(eligible) - 1L)
  }
  test <- with_seed(seed, eligible[sample.int(length(eligible), n_test)])
  train <- setdiff(perts, test)

  if (type == "unseen_single") {
    train_parts <- combo_parts(train, x$combo_sep)
    leaks <- vapply(train_parts, function(g) any(g %in% test), logical(1))
    train <- train[!leaks]
  }

  structure(
    list(train = train, test = test, control = x$control, type = type,
         test_frac = test_frac, seed = seed),
    class = "pe_split"
  )
}

#' @export
print.pe_split <- function(x, ...) {
  cat("<pe_split: ", x$type, ">\n", sep = "")
  cat("  train: ", length(x$train), " perturbations\n", sep = "")
  cat("  test:  ", length(x$test), " perturbations (",
      paste(utils::head(x$test, 5), collapse = ", "),
      if (length(x$test) > 5) ", ..." else "", ")\n", sep = "")
  invisible(x)
}
