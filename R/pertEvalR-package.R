#' pertEvalR: Calibrated Evaluation of Single-Cell Perturbation Predictions
#'
#' Tools for evaluating models that predict transcriptional responses to
#' perturbations. The typical workflow is:
#'
#' 1. Build a dataset with [pe_data()] or [as_pe_data()].
#' 2. Define a held-out set of perturbations with [pe_split()].
#' 3. Produce predictions, either from your own model or from a reference
#'    baseline such as [baseline_control()] or [baseline_mean()].
#' 4. Compare predictions against [pe_truth()] (metrics are planned; see the
#'    package README for the roadmap).
#'
#' @keywords internal
"_PACKAGE"
