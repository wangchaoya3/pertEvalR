# These tests run only when the optional packages are installed (always on the
# GitHub Actions check, which installs all Suggests).

small_sim <- function() {
  pe_simulate(n_genes = 30, n_single = 3, n_combo = 1, cells_per_pert = 5,
              n_control = 10, n_affected = 5, seed = 1)
}

test_that("as_pe_data works for SingleCellExperiment objects", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  x <- small_sim()
  cd <- data.frame(target = x$perturbation, row.names = colnames(x$expr))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = x$expr), colData = cd
  )

  y <- as_pe_data(sce, perturbation_col = "target")
  expect_s3_class(y, "pe_data")
  expect_identical(y$perturbation, x$perturbation)
  expect_equal(pe_pseudobulk(y), pe_pseudobulk(x))
  expect_error(as_pe_data(sce, perturbation_col = "missing"), "not found")
})

test_that("as_pe_data works for Seurat objects", {
  skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
  x <- small_sim()
  # Seurat expects non-negative values in the counts layer.
  counts <- Matrix::Matrix(x$expr - min(x$expr), sparse = TRUE)
  seu <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  seu <- SeuratObject::AddMetaData(
    seu, metadata = stats::setNames(x$perturbation, colnames(counts)),
    col.name = "target"
  )

  y <- as_pe_data(seu, perturbation_col = "target", layer = "counts")
  ref <- pe_data(counts, x$perturbation)
  expect_s3_class(y, "pe_data")
  expect_identical(y$perturbation, x$perturbation)
  expect_equal(pe_pseudobulk(y), pe_pseudobulk(ref))
  expect_error(as_pe_data(seu, perturbation_col = "missing"), "not found")
})
