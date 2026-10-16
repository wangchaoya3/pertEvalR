#' Wrap predictions for evaluation
#'
#' Turns a genes x perturbations matrix of predicted pseudobulk expression,
#' for example the output of your own model, into a `pe_prediction` object
#' that [pe_evaluate()] understands. The reference baselines return objects of
#' the same class.
#'
#' Predictions must be on the same scale as the data in the `pe_data` object
#' (for example log-normalized expression averaged over cells).
#'
#' @param pred A numeric matrix (or data frame, or \pkg{Matrix} matrix) with
#'   gene names as row names and one column per held-out perturbation of
#'   `split`, named by perturbation label. Columns may be in any order; they
#'   are reordered to match `split$test`. Extra genes are allowed and ignored
#'   during evaluation.
#' @param split The `pe_split` object the predictions were made for.
#' @param method A short name for the model, used to label results.
#'
#' @return An object of class `pe_prediction`: a list with `pred`, `method`
#'   and `split`.
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' s <- pe_split(x, seed = 1)
#' # A (cheating) "oracle" model that returns the observed profiles
#' oracle <- pe_prediction(pe_truth(x, s), s, method = "oracle")
#' oracle
pe_prediction <- function(pred, split, method = "custom") {
  if (!inherits(split, "pe_split")) {
    stop("`split` must be a `pe_split` object; see ?pe_split.", call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("`method` must be a single string.", call. = FALSE)
  }
  pred <- as.matrix(pred)
  if (!is.numeric(pred)) {
    stop("`pred` must be numeric.", call. = FALSE)
  }
  if (is.null(rownames(pred)) || anyDuplicated(rownames(pred)) > 0) {
    stop("`pred` must have unique gene names as row names.", call. = FALSE)
  }
  cols <- colnames(pred)
  if (is.null(cols) || anyDuplicated(cols) > 0 || !setequal(cols, split$test)) {
    stop("The columns of `pred` must be named after exactly the held-out ",
         "perturbations of `split` (split$test).", call. = FALSE)
  }
  pred <- pred[, split$test, drop = FALSE]
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
