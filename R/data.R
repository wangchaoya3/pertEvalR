#' Create a perturbation dataset
#'
#' `pe_data()` bundles a genes x cells expression matrix with the
#' perturbation label of every cell. It is the input to every other function
#' in the package.
#'
#' @param expr A genes x cells expression matrix, either a base `matrix` or a
#'   sparse matrix from the \pkg{Matrix} package. Values should already be
#'   normalized (for example log-normalized counts); `pe_data()` does not
#'   transform them. Row names (gene identifiers) are required.
#' @param perturbation A character or factor vector with one label per cell
#'   (per column of `expr`).
#' @param control The label used for unperturbed control cells.
#' @param combo_sep The separator used in labels of combinatorial
#'   perturbations, for example `"+"` in `"GENE1+GENE2"`.
#'
#' @return An object of class `pe_data`: a list with elements `expr`,
#'   `perturbation`, `control` and `combo_sep`.
#' @seealso [as_pe_data()] to build one from a `SingleCellExperiment` or
#'   `Seurat` object; [pe_simulate()] for a toy dataset.
#' @export
#' @examples
#' expr <- matrix(rnorm(40), nrow = 4,
#'                dimnames = list(paste0("g", 1:4), NULL))
#' pert <- rep(c("control", "A", "B", "A+B", "control"), each = 2)
#' x <- pe_data(expr, pert)
#' x
pe_data <- function(expr, perturbation, control = "control", combo_sep = "+") {
  if (!is_matrix_like(expr)) {
    stop("`expr` must be a matrix or a Matrix-package sparse matrix.",
         call. = FALSE)
  }
  if (is.null(rownames(expr)) || anyDuplicated(rownames(expr)) > 0) {
    stop("`expr` must have unique gene names as row names.", call. = FALSE)
  }
  if (!is.character(control) || length(control) != 1 || is.na(control)) {
    stop("`control` must be a single label.", call. = FALSE)
  }
  if (!is.character(combo_sep) || length(combo_sep) != 1 || !nzchar(combo_sep)) {
    stop("`combo_sep` must be a single non-empty string.", call. = FALSE)
  }
  perturbation <- as.character(perturbation)
  if (length(perturbation) != ncol(expr)) {
    stop("`perturbation` has ", length(perturbation), " entries but `expr` has ",
         ncol(expr), " cells (columns).", call. = FALSE)
  }
  if (anyNA(perturbation)) {
    stop("`perturbation` must not contain missing values.", call. = FALSE)
  }
  if (!control %in% perturbation) {
    stop("No cells are labelled with the control label \"", control, "\".",
         call. = FALSE)
  }
  structure(
    list(expr = expr, perturbation = perturbation,
         control = control, combo_sep = combo_sep),
    class = "pe_data"
  )
}

#' Convert single-cell objects to a perturbation dataset
#'
#' Methods are provided for `SingleCellExperiment` and `Seurat` (v5) objects.
#' The default method passes its arguments on to [pe_data()].
#'
#' @param x A `SingleCellExperiment`, a `Seurat` object, or (for the default
#'   method) an expression matrix.
#' @param ... Passed on to methods. For the default method, the remaining
#'   arguments of [pe_data()].
#' @param perturbation_col Name of the cell-metadata column (`colData` or
#'   `meta.data`) holding the perturbation labels.
#' @param control,combo_sep See [pe_data()].
#' @param assay For `SingleCellExperiment`: the assay to use
#'   (default `"logcounts"`). For `Seurat`: the assay name (default: the
#'   object's default assay).
#' @param layer For `Seurat`: the layer to use (default `"data"`).
#'
#' @return A `pe_data` object.
#' @export
as_pe_data <- function(x, ...) {
  UseMethod("as_pe_data")
}

#' @rdname as_pe_data
#' @export
as_pe_data.default <- function(x, ...) {
  pe_data(x, ...)
}

#' @rdname as_pe_data
#' @export
as_pe_data.SingleCellExperiment <- function(x, perturbation_col,
                                            control = "control",
                                            combo_sep = "+",
                                            assay = "logcounts", ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.", call. = FALSE)
  }
  meta <- SummarizedExperiment::colData(x)
  if (!perturbation_col %in% colnames(meta)) {
    stop("Column \"", perturbation_col, "\" not found in colData(x).",
         call. = FALSE)
  }
  pe_data(SummarizedExperiment::assay(x, assay),
          perturbation = meta[[perturbation_col]],
          control = control, combo_sep = combo_sep)
}

#' @rdname as_pe_data
#' @export
as_pe_data.Seurat <- function(x, perturbation_col,
                              control = "control",
                              combo_sep = "+",
                              assay = NULL, layer = "data", ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required.", call. = FALSE)
  }
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(x)
  }
  meta <- x[[]]
  if (!perturbation_col %in% colnames(meta)) {
    stop("Column \"", perturbation_col, "\" not found in the Seurat meta.data.",
         call. = FALSE)
  }
  expr <- SeuratObject::LayerData(x, assay = assay, layer = layer)
  pe_data(expr, perturbation = meta[[perturbation_col]],
          control = control, combo_sep = combo_sep)
}

#' List the perturbations in a dataset
#'
#' @param x A `pe_data` object.
#' @param include_control Whether to include the control label.
#' @param type Which perturbations to return: `"all"`, `"single"` or `"combo"`.
#'
#' @return A character vector of perturbation labels, in order of first
#'   appearance.
#' @export
#' @examples
#' x <- pe_simulate(seed = 1)
#' head(pe_perturbations(x, type = "combo"))
pe_perturbations <- function(x, include_control = FALSE,
                             type = c("all", "single", "combo")) {
  check_pe_data(x)
  type <- match.arg(type)
  perts <- unique(x$perturbation)
  ctrl <- perts == x$control
  combo <- is_combo(perts, x$combo_sep) & !ctrl
  keep <- switch(type,
    all = !ctrl,
    single = !ctrl & !combo,
    combo = combo
  )
  if (include_control) {
    keep <- keep | ctrl
  }
  perts[keep]
}

#' @export
print.pe_data <- function(x, ...) {
  n_single <- length(pe_perturbations(x, type = "single"))
  n_combo <- length(pe_perturbations(x, type = "combo"))
  n_ctrl <- sum(x$perturbation == x$control)
  cat("<pe_data>\n")
  cat("  genes:         ", nrow(x$expr), "\n", sep = "")
  cat("  cells:         ", ncol(x$expr), " (", n_ctrl, " control)\n", sep = "")
  cat("  perturbations: ", n_single, " single, ", n_combo, " combinatorial\n",
      sep = "")
  cat("  control label: \"", x$control, "\"\n", sep = "")
  invisible(x)
}
