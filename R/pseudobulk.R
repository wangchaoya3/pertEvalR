#' Average expression per perturbation (pseudobulk)
#'
#' Computes the mean expression profile of the cells in each perturbation
#' group. This is the representation used by most perturbation-prediction
#' metrics and by the reference baselines.
#'
#' @param x A `pe_data` object.
#' @param perturbations Labels of the groups to aggregate. Defaults to all
#'   labels in `x`, including the control.
#'
#' @return A dense genes x perturbations numeric matrix, with columns in the
#'   order of `perturbations`.
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' pb <- pe_pseudobulk(x, c("control", "P1"))
#' dim(pb)
pe_pseudobulk <- function(x, perturbations = NULL) {
  check_pe_data(x)
  if (is.null(perturbations)) {
    perturbations <- unique(x$perturbation)
  }
  perturbations <- as.character(perturbations)
  if (length(perturbations) == 0) {
    stop("`perturbations` must contain at least one label.", call. = FALSE)
  }
  if (anyDuplicated(perturbations) > 0) {
    stop("`perturbations` must not contain duplicates.", call. = FALSE)
  }
  missing <- setdiff(perturbations, x$perturbation)
  if (length(missing) > 0) {
    stop("Unknown perturbation(s): ",
         paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
  }

  # Sparse cells x groups averaging matrix: column j holds 1/n_j for the n_j
  # cells of group j, so expr %*% avg gives the group means.
  group <- match(x$perturbation, perturbations)
  cells <- which(!is.na(group))
  n_cells <- tabulate(group[cells], nbins = length(perturbations))
  avg <- Matrix::sparseMatrix(
    i = cells,
    j = group[cells],
    x = 1 / n_cells[group[cells]],
    dims = c(ncol(x$expr), length(perturbations))
  )
  out <- as.matrix(x$expr %*% avg)
  dimnames(out) <- list(rownames(x$expr), perturbations)
  out
}

#' Observed pseudobulk profiles of the held-out perturbations
#'
#' Convenience wrapper returning the ground truth that predictions for a split
#' are compared against.
#'
#' @param x A `pe_data` object.
#' @param split A `pe_split` object created from `x`.
#'
#' @return A genes x test-perturbations matrix (see [pe_pseudobulk()]).
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' dim(pe_truth(x, s))
pe_truth <- function(x, split) {
  check_pe_data(x)
  check_pe_split(split, x)
  pe_pseudobulk(x, split$test)
}
