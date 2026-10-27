# Validation against an existing tool, and runtime
# ---------------------------------------------------------------------------
# Part A. Reference-implementation check (Norman et al. 2019, unseen singles)
#   pertEvalR's baseline_linear() re-implements the linear baseline of
#   Ahlmann-Eltze, Huber and Anders (2025, Nature Methods). Here it is compared
#   with the authors' own code (repository const-ae/linear_perturbation_
#   prediction-Paper, MIT licence, file benchmark/src/run_linear_pretrained_model.R):
#   their function solve_y_axb() is read from that file and run with their
#   embedding and training steps (PCA of the training pseudobulk with irlba,
#   control included in training, gene embedding = PCA scores, perturbation
#   embedding = the target gene's row, zero for control; pca_dim = 10,
#   ridge = 0.1). Both implementations receive the same pseudobulk input and
#   the same five splits, so differences come from the model code only.
#   Reported: per-perturbation correlation between the two predicted changes
#   from control, agreement of the per-perturbation scores, paired differences
#   (pe_compare), and whether both reproduce the published finding that the
#   linear model is close to the training mean.
#
# Part B. Runtime and memory of the full workflow on Norman (all perturbations)
#   and on Frangieh et al. 2021 (IFN-gamma condition), one split each.
#
# Needs: the Norman cache (real02), the Frangieh cache (real05; Part B only,
# skipped if absent) and the irlba package (install.packages("irlba")).
# The reference file is downloaded once to analysis/data/; if the download
# fails (network), save the file there by hand from
#   https://github.com/const-ae/linear_perturbation_prediction-Paper/blob/main/benchmark/src/run_linear_pretrained_model.R
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real09_validation.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
devtools::load_all(pkg, quiet = TRUE)
if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Please install irlba first: install.packages(\"irlba\")")
}

data_dir <- file.path(pkg, "analysis", "data")
out_dir <- file.path(pkg, "analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seeds <- 1:5
pca_dim <- 10
ridge <- 0.1

# A1. The authors' solve_y_axb() ----------------------------------------------
ref_file <- file.path(data_dir, "run_linear_pretrained_model.R")
ref_url <- paste0("https://raw.githubusercontent.com/const-ae/",
                  "linear_perturbation_prediction-Paper/main/benchmark/src/",
                  "run_linear_pretrained_model.R")
if (!file.exists(ref_file)) {
  ok <- tryCatch({
    utils::download.file(ref_url, ref_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(ref_file)) {
    stop("Could not download the reference file. Save it by hand as\n  ",
         ref_file, "\nfrom ", ref_url)
  }
}
# Evaluate only the definition of solve_y_axb (the rest of the file is a
# command-line script that would try to read its own arguments and data).
exprs <- parse(ref_file, keep.source = FALSE)
is_def <- vapply(exprs, function(e) {
  is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
    identical(as.character(e[[2]]), "solve_y_axb")
}, logical(1))
if (sum(is_def) != 1) stop("solve_y_axb() not found in ", ref_file)
ref_env <- new.env(parent = globalenv())
eval(exprs[[which(is_def)]], envir = ref_env)
solve_y_axb <- ref_env$solve_y_axb
cat("Reference solve_y_axb() loaded; arguments:",
    paste(names(formals(solve_y_axb)), collapse = ", "), "\n")

# The authors' prediction steps, applied to pertEvalR's pseudobulk
reference_linear <- function(x, split) {
  train <- split$train[!is_combo(split$train, x$combo_sep)]
  test <- split$test
  ctrl_label <- x$control
  X <- pe_pseudobulk(x, c(train, ctrl_label))
  baseline <- X[, ctrl_label]
  Y <- X - baseline
  pca <- irlba::prcomp_irlba(X, n = pca_dim)
  gene_emb <- pca$x
  rownames(gene_emb) <- rownames(X)
  pert_emb <- t(gene_emb)                      # columns = genes
  emb_of <- function(labels) {
    m <- matrix(0, nrow = pca_dim, ncol = length(labels),
                dimnames = list(NULL, labels))
    real <- labels != ctrl_label
    m[, real] <- pert_emb[, labels[real], drop = FALSE]
    m
  }
  coefs <- solve_y_axb(Y = Y, A = gene_emb, B = emb_of(colnames(X)),
                       A_ridge = ridge, B_ridge = ridge)
  pred <- gene_emb %*% coefs$K %*% emb_of(test) + coefs$center + baseline
  pred <- as.matrix(pred)
  dimnames(pred) <- list(rownames(X), test)
  pe_prediction(pred, split, method = "reference_linear")
}

# A2. Norman, five splits ----------------------------------------------------------
norman_cache <- file.path(data_dir, "norman_pe_data.rds")
if (!file.exists(norman_cache)) stop("Run real02_norman_baselines.R first.")
x <- readRDS(norman_cache)
ctrl <- pe_pseudobulk(x, x$control)[, 1]

agree <- list()
scores <- list()
for (sd in seeds) {
  s <- pe_split(x, type = "unseen_single", test_frac = 0.25, seed = sd)
  p_pe <- baseline_linear(x, s, n_comp = pca_dim, ridge = ridge)
  p_ref <- reference_linear(x, s)
  p_mean <- baseline_mean(x, s)
  a <- p_pe$pred - ctrl
  b <- p_ref$pred[rownames(a), colnames(a)] - ctrl
  agree[[sd]] <- data.frame(
    seed = sd, perturbation = colnames(a),
    cor_change = vapply(seq_len(ncol(a)), function(j) cor(a[, j], b[, j]),
                        numeric(1)),
    rel_diff = sqrt(colSums((a - b)^2) / colSums(b^2)))
  r <- pe_evaluate(list(p_mean, p_pe, p_ref), x)
  r$seed <- sd
  scores[[sd]] <- r
  cat(sprintf("Split %d/%d done\n", sd, length(seeds)))
}
agree <- do.call(rbind, agree)
res <- do.call(rbind, scores)
res$perturbation <- paste0("split", res$seed, ":", res$perturbation)
write.csv(agree, file.path(out_dir, "validation_reference_agreement.csv"),
          row.names = FALSE)
write.csv(res, file.path(out_dir, "validation_reference_scores.csv"),
          row.names = FALSE)

cat("\n== Predicted changes from control: pertEvalR vs reference code ==\n")
cat(sprintf("held-out perturbations (5 splits): %d\n", nrow(agree)))
print(round(quantile(agree$cor_change, c(0, 0.05, 0.5, 1)), 4))
cat("relative difference ||a - b|| / ||b||, median and max:",
    signif(median(agree$rel_diff), 3), signif(max(agree$rel_diff), 3), "\n")

# Score agreement per metric
w <- reshape(res[res$method %in% c("linear", "reference_linear"),
                 c("method", "perturbation", "metric", "value")],
             idvar = c("perturbation", "metric"), timevar = "method",
             direction = "wide")
score_agree <- do.call(rbind, lapply(split(w, w$metric), function(d) {
  ok <- stats::complete.cases(d)
  data.frame(metric = d$metric[1], n = sum(ok),
             cor = cor(d$value.linear[ok], d$value.reference_linear[ok]),
             mean_pertEvalR = mean(d$value.linear[ok]),
             mean_reference = mean(d$value.reference_linear[ok]),
             max_abs_diff = max(abs(d$value.linear[ok] -
                                      d$value.reference_linear[ok])))
}))
rownames(score_agree) <- NULL
write.csv(score_agree, file.path(out_dir, "validation_score_agreement.csv"),
          row.names = FALSE)
cat("\n== Per-perturbation scores: pertEvalR linear vs reference code ==\n")
print(transform(score_agree, cor = round(cor, 4),
                mean_pertEvalR = signif(mean_pertEvalR, 4),
                mean_reference = signif(mean_reference, 4),
                max_abs_diff = signif(max_abs_diff, 3)), row.names = FALSE)

cmp <- pe_compare(res[c("method", "perturbation", "metric", "value")],
                  reference = "mean", seed = 1)
write.csv(cmp, file.path(out_dir, "validation_reference_compare.csv"),
          row.names = FALSE)
cat("\n== Each linear implementation vs training mean (diff = method - mean) ==\n")
show <- cmp[, c("metric", "method_a", "n", "diff", "ci_low", "ci_high",
                "favours")]
show[c("diff", "ci_low", "ci_high")] <- round(show[c("diff", "ci_low",
                                                     "ci_high")], 4)
print(show[order(show$metric, show$method_a), ], row.names = FALSE)

# B. Runtime and memory --------------------------------------------------------------
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, sec = proc.time()[["elapsed"]] - t0)
}
run_workflow <- function(x, label) {
  st <- list()
  add <- function(step, tm) {
    st[[length(st) + 1]] <<- data.frame(dataset = label, step = step,
                                        seconds = round(tm, 2))
  }
  r <- timed(pe_signal(x, seed = 1)); add("pe_signal", r$sec)
  s <- pe_split(x, type = "unseen_single", test_frac = 0.25, seed = 1)
  r <- timed(list(baseline_control(x, s), baseline_mean(x, s),
                  baseline_linear(x, s)))
  add("baselines (no change, mean, linear)", r$sec)
  preds <- r$value
  r <- timed(pe_evaluate(preds, x)); add("pe_evaluate", r$sec)
  ev <- r$value
  r <- timed(pe_bounds(x, s, seed = 1)); add("pe_bounds", r$sec)
  bd <- r$value
  r <- timed(pe_calibrate(ev, bd)); add("pe_calibrate", r$sec)
  r <- timed(pe_compare(ev, reference = "mean", seed = 1))
  add("pe_compare (2000 bootstrap)", r$sec)
  out <- do.call(rbind, st)
  n_single <- sum(!is_combo(pe_perturbations(x), x$combo_sep))
  attr(out, "size") <- data.frame(
    dataset = label, genes = nrow(x$expr), cells = ncol(x$expr),
    perturbations = length(pe_perturbations(x)),
    single = n_single, test_perturbations = length(s$test),
    object_MB = round(as.numeric(utils::object.size(x)) / 2^20, 1),
    total_seconds = round(sum(out$seconds), 1))
  out
}

rt <- list(run_workflow(x, "Norman"))
fr_cache <- file.path(data_dir, "frangieh_cache.rds")
if (file.exists(fr_cache)) {
  fr <- readRDS(fr_cache)
  cells <- fr$condition == "IFNg"
  pert <- fr$perturbation[cells]
  n <- table(pert)
  ok <- names(n)[n >= 30 | names(n) == "control"]
  keep <- which(cells)[pert %in% ok]
  xf <- pe_data(fr$expr[, keep], fr$perturbation[keep], control = "control")
  rm(fr); invisible(gc(verbose = FALSE))
  rt[[2]] <- run_workflow(xf, "Frangieh IFN-gamma")
} else {
  cat("Frangieh cache not found; runtime only for Norman.\n")
}
runtime <- do.call(rbind, rt)
sizes <- do.call(rbind, lapply(rt, attr, "size"))
write.csv(runtime, file.path(out_dir, "validation_runtime.csv"),
          row.names = FALSE)
write.csv(sizes, file.path(out_dir, "validation_runtime_sizes.csv"),
          row.names = FALSE)
cat("\n== Runtime (seconds, one split) ==\n")
print(runtime, row.names = FALSE)
print(sizes, row.names = FALSE)

cat("\n== Machine ==\n")
cat("R:", R.version.string, "\n")
cat("OS:", utils::sessionInfo()$running, "\n")
cat("CPU:", Sys.getenv("PROCESSOR_IDENTIFIER"), "| logical cores:",
    parallel::detectCores(), "\n")
if (.Platform$OS.type == "windows") {
  mem <- tryCatch(system("wmic ComputerSystem get TotalPhysicalMemory",
                         intern = TRUE), error = function(e) "")
  mem <- suppressWarnings(as.numeric(gsub("\\D", "", mem)))
  mem <- mem[!is.na(mem)]
  if (length(mem)) cat("RAM:", round(mem[1] / 2^30, 1), "GB\n")
}
cat("pertEvalR", as.character(utils::packageVersion("pertEvalR")),
    "| irlba", as.character(utils::packageVersion("irlba")), "\n")
