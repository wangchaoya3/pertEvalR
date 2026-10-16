toy_expr <- function() {
  matrix(as.numeric(1:20), nrow = 2,
         dimnames = list(c("g1", "g2"), paste0("c", 1:10)))
}
toy_pert <- c(rep("control", 4), "A", "A", "B", "B", "A+B", "A+B")

test_that("pe_data validates its inputs", {
  expr <- toy_expr()
  expect_s3_class(pe_data(expr, toy_pert), "pe_data")

  expect_error(pe_data(as.data.frame(expr), toy_pert), "matrix")
  expect_error(pe_data(unname(expr), toy_pert), "row names")
  expect_error(pe_data(expr, toy_pert[-1]), "entries")
  expect_error(pe_data(expr, replace(toy_pert, 1, NA)), "missing")
  expect_error(pe_data(expr, toy_pert, control = "ctrl"), "control label")
})

test_that("pe_data accepts sparse matrices", {
  sparse <- Matrix::Matrix(toy_expr(), sparse = TRUE)
  expect_s3_class(pe_data(sparse, toy_pert), "pe_data")
})

test_that("as_pe_data default method forwards to pe_data", {
  x <- as_pe_data(toy_expr(), perturbation = toy_pert)
  expect_identical(x, pe_data(toy_expr(), toy_pert))
})

test_that("pe_perturbations separates singles and combinations", {
  x <- pe_data(toy_expr(), toy_pert)
  expect_identical(pe_perturbations(x), c("A", "B", "A+B"))
  expect_identical(pe_perturbations(x, type = "single"), c("A", "B"))
  expect_identical(pe_perturbations(x, type = "combo"), "A+B")
  expect_identical(pe_perturbations(x, include_control = TRUE),
                   c("control", "A", "B", "A+B"))
})

test_that("pe_pseudobulk returns per-group means", {
  x <- pe_data(toy_expr(), toy_pert)
  pb <- pe_pseudobulk(x, c("A", "control"))
  expect_identical(dimnames(pb), list(c("g1", "g2"), c("A", "control")))
  expected_a <- rowMeans(toy_expr()[, toy_pert == "A"])
  expected_ctrl <- rowMeans(toy_expr()[, toy_pert == "control"])
  expect_equal(pb[, "A"], expected_a)
  expect_equal(pb[, "control"], expected_ctrl)
})

test_that("pe_pseudobulk gives the same result for dense and sparse input", {
  dense <- pe_data(toy_expr(), toy_pert)
  sparse <- pe_data(Matrix::Matrix(toy_expr(), sparse = TRUE), toy_pert)
  expect_equal(pe_pseudobulk(dense), pe_pseudobulk(sparse))
})

test_that("pe_pseudobulk rejects unknown or duplicated labels", {
  x <- pe_data(toy_expr(), toy_pert)
  expect_error(pe_pseudobulk(x, "Z"), "Unknown")
  expect_error(pe_pseudobulk(x, c("A", "A")), "duplicates")
})

test_that("print.pe_data reports the dataset shape", {
  x <- pe_data(toy_expr(), toy_pert)
  expect_output(print(x), "2 single, 1 combinatorial")
})
