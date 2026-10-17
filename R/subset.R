#' Restrict results or bounds to a set of perturbations
#'
#' Keeps only the rows of a [pe_evaluate()] result or a [pe_bounds()] table
#' that belong to the given perturbations, while preserving the evaluation
#' settings that [pe_calibrate()] checks. Plain row subsetting with `[` drops
#' those settings.
#'
#' A typical use is to evaluate only perturbations with a measurable effect,
#' for example those whose split-half reliability from [pe_signal()] exceeds
#' a threshold chosen before looking at any model. On datasets where most
#' perturbations barely change expression, averaging over all of them mainly
#' measures how models fit noise.
#'
#' Metrics that compare perturbations with each other (`pds`) keep the values
#' computed against all held-out perturbations; subsetting only changes which
#' perturbations are averaged.
#'
#' @param x A data frame returned by [pe_evaluate()] or [pe_bounds()].
#' @param perturbations Character vector of perturbations to keep, or a
#'   `pe_signal` object, in which case the perturbations with reliability of
#'   at least `min_reliability` are kept.
#' @param min_reliability Reliability threshold used when `perturbations` is a
#'   `pe_signal` object.
#'
#' @return `x` restricted to the selected perturbations, with its attributes
#'   kept. Selected perturbations absent from `x` are ignored; it is an error
#'   if none of them is present.
#' @seealso [pe_signal()], [pe_calibrate()]
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' res <- pe_evaluate(baseline_mean(x, s), x)
#' sig <- pe_signal(x, seed = 1)
#' res_rel <- pe_subset(res, sig, min_reliability = 0.5)
#' b_rel <- pe_subset(pe_bounds(x, s, seed = 1), sig, min_reliability = 0.5)
#' pe_calibrate(res_rel, b_rel)
pe_subset <- function(x, perturbations, min_reliability = 0.5) {
  if (!is.data.frame(x) || !"perturbation" %in% names(x)) {
    stop("`x` must be a data frame returned by pe_evaluate() or pe_bounds().",
         call. = FALSE)
  }
  if (inherits(perturbations, "pe_signal")) {
    if (!is.numeric(min_reliability) || length(min_reliability) != 1 ||
        is.na(min_reliability)) {
      stop("`min_reliability` must be a single number.", call. = FALSE)
    }
    per <- perturbations$per_perturbation
    perturbations <- per$perturbation[!is.na(per$reliability) &
                                        per$reliability >= min_reliability]
  }
  if (!is.character(perturbations)) {
    stop("`perturbations` must be a character vector or a pe_signal object.",
         call. = FALSE)
  }
  keep <- x$perturbation %in% perturbations
  if (!any(keep)) {
    stop("None of the selected perturbations is present in `x`.",
         call. = FALSE)
  }
  attrs <- attributes(x)
  out <- x[keep, , drop = FALSE]
  rownames(out) <- NULL
  extra <- setdiff(names(attrs), c("names", "row.names", "class"))
  for (a in extra) {
    attr(out, a) <- attrs[[a]]
  }
  out
}
