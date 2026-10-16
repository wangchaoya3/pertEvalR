# Internal helpers --------------------------------------------------------

# Evaluate `expr` with a temporary RNG seed, restoring the caller's RNG state
# afterwards. `expr` is evaluated lazily, after `set.seed()` has run.
with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

# TRUE for labels describing a combinatorial perturbation ("A+B").
is_combo <- function(labels, sep) {
  grepl(sep, labels, fixed = TRUE)
}

# Split perturbation labels into their component genes.
combo_parts <- function(labels, sep) {
  strsplit(labels, sep, fixed = TRUE)
}

check_pe_data <- function(x, arg = "x") {
  if (!inherits(x, "pe_data")) {
    stop("`", arg, "` must be a `pe_data` object; see ?pe_data.", call. = FALSE)
  }
  invisible(x)
}

check_pe_split <- function(split, x) {
  if (!inherits(split, "pe_split")) {
    stop("`split` must be a `pe_split` object; see ?pe_split.", call. = FALSE)
  }
  if (!identical(split$control, x$control)) {
    stop("`split` was created for a dataset with a different control label.",
         call. = FALSE)
  }
  unknown <- setdiff(c(split$train, split$test), x$perturbation)
  if (length(unknown) > 0) {
    stop("`split` refers to perturbations not present in `x`: ",
         paste(utils::head(unknown, 5), collapse = ", "), call. = FALSE)
  }
  invisible(split)
}

is_matrix_like <- function(expr) {
  is.matrix(expr) || methods::is(expr, "Matrix")
}
