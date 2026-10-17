# Frangieh et al. 2021: baselines, calibration and signal per culture condition
# ---------------------------------------------------------------------------
# Data: scPerturb FrangiehIzar2021_RNA.h5ad (Zenodo 10044268, CC BY 4.0).
# Melanoma cells, CRISPR knockouts of ~250 genes, profiled under three culture
# conditions (obs column `perturbation_2`): monoculture ("Control"), IFN-gamma
# stimulation ("IFNg") and co-culture with T cells ("Co-culture").
#
# Cell filtering (decided from the metadata inspection):
#   * only cells with exactly one detected guide (MOI == 1). Cells with MOI 0
#     carry no detected guide but are labelled "control" in scPerturb; cells
#     with MOI >= 2 carry several guides but only one perturbation label.
#     Neither is a clean single-knockout or non-targeting cell.
#   * perturbations whose target gene is not measured are dropped.
#   * each condition is analysed separately, against its own control cells;
#     perturbations with fewer than `min_cells` cells in a condition are
#     dropped from that condition.
#
# Per condition: pe_signal(); 5 random unseen_single splits with the no-change,
# training-mean and linear baselines; pe_compare() vs the training mean;
# pe_bounds() + pe_calibrate().
#
# Memory: the full matrix is ~8 GB, but it is read in gene blocks and only the
# kept cells and selected genes are stored (~2-3 GB). The first run builds a
# cache (analysis/data/frangieh_cache.rds); later runs start from it.
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real05_frangieh_baselines.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
devtools::load_all(pkg, quiet = TRUE)

data_dir <- file.path(pkg, "analysis", "data")
out_dir <- file.path(pkg, "analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
h5ad <- file.path(data_dir, "FrangiehIzar2021_RNA.h5ad")
cache <- file.path(data_dir, "frangieh_cache.rds")

# Settings
n_hvg <- 5000          # highly variable genes kept (plus all target genes)
max_cells_pert <- 300  # cells kept per perturbation and condition
max_cells_ctrl <- 3000 # control cells kept per condition
min_cells <- 30        # minimum cells for a perturbation within a condition
block_genes <- 1500    # genes read from disk at a time
seeds <- 1:5           # random splits per condition

# Short ASCII names for the three conditions
condition_name <- function(v) {
  v <- as.character(v)
  out <- rep(NA_character_, length(v))
  out[v == "Control"] <- "monoculture"
  out[grepl("IFN", v, fixed = TRUE)] <- "IFNg"
  out[v == "Co-culture"] <- "coculture"
  if (anyNA(out)) stop("Unexpected condition values: ",
                       paste(unique(v[is.na(out)]), collapse = ", "))
  out
}

# 1. Build (or load) the expression cache ------------------------------------------

build_frangieh <- function() {
  enc <- rhdf5::h5readAttributes(h5ad, "/X")[["encoding-type"]]
  cat("X encoding:", enc, "\n")
  if (!identical(enc, "csc_matrix")) {
    stop("Expected X stored as csc_matrix (by gene); found ", enc, ".")
  }

  ad <- anndataR::read_h5ad(h5ad, as = "HDF5AnnData")
  obs <- ad$obs
  genes <- make.unique(rownames(ad$var))
  rm(ad)
  rhdf5::h5closeAll()

  pert <- as.character(obs$perturbation)
  cond <- condition_name(obs$perturbation_2)
  libsize <- obs$ncounts
  n_cells <- length(pert)

  # Clean cells: exactly one detected guide
  clean <- obs$MOI == 1
  cat("Cells with exactly one guide:", sum(clean), "of", n_cells, "\n")
  ctrl_guides <- as.character(obs$guide_id)[clean & pert == "control"]
  cat("Control guide families:\n")
  print(table(sub("_[0-9]+$", "", ctrl_guides)))

  # Drop perturbations whose target gene is not measured
  labels <- setdiff(unique(pert[clean]), "control")
  absent <- labels[!labels %in% genes]
  if (length(absent) > 0) {
    cat("Dropping", length(absent), "perturbation(s) with unmeasured targets:",
        paste(absent, collapse = ", "), "\n")
  }
  clean <- clean & !pert %in% absent
  targets <- setdiff(labels, absent)

  # Cap cells per perturbation within each condition
  set.seed(1)
  idx_clean <- which(clean)
  groups <- split(idx_clean, paste(cond[idx_clean], pert[idx_clean], sep = "|"))
  keep <- sort(unlist(lapply(names(groups), function(k) {
    idx <- groups[[k]]
    cap <- if (endsWith(k, "|control")) max_cells_ctrl else max_cells_pert
    if (length(idx) > cap) idx[sample.int(length(idx), cap)] else idx
  }), use.names = FALSE))
  n_keep <- length(keep)
  cell_map <- integer(n_cells)
  cell_map[keep] <- seq_len(n_keep)
  scale <- 1e4 / libsize
  cat("Cells kept:", n_keep, "\n")
  print(table(condition = cond[keep], control = pert[keep] == "control"))

  indptr <- as.numeric(rhdf5::h5read(h5ad, "/X/indptr"))
  stopifnot(length(indptr) == length(genes) + 1)

  # Read genes g_from..g_to; return log1p(CP10k) values of kept cells
  read_block <- function(g_from, g_to) {
    s <- indptr[g_from]
    n <- indptr[g_to + 1] - s
    if (n == 0) return(NULL)
    v <- as.numeric(rhdf5::h5read(h5ad, "/X/data", start = s + 1, count = n))
    i <- as.integer(rhdf5::h5read(h5ad, "/X/indices", start = s + 1,
                                  count = n)) + 1L
    g <- rep.int(g_from:g_to, diff(indptr[g_from:(g_to + 1)]))
    ok <- cell_map[i] > 0
    list(cell = cell_map[i[ok]], gene = g[ok],
         value = log1p(v[ok] * scale[i[ok]]))
  }

  blocks <- split(seq_along(genes), ceiling(seq_along(genes) / block_genes))

  # Pass 1: per-gene mean and variance over kept cells
  s1 <- numeric(length(genes))
  s2 <- numeric(length(genes))
  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    r <- read_block(min(b), max(b))
    if (!is.null(r) && length(r$value) > 0) {
      agg <- rowsum(cbind(r$value, r$value^2), r$gene)
      idx <- as.integer(rownames(agg))
      s1[idx] <- agg[, 1]
      s2[idx] <- agg[, 2]
    }
    cat(sprintf("\rPass 1 (gene variance): block %d/%d", k, length(blocks)))
    rm(r); invisible(gc(verbose = FALSE))
  }
  cat("\n")
  mu <- s1 / n_keep
  vr <- s2 / n_keep - mu^2
  hvg <- order(vr, decreasing = TRUE)[seq_len(min(n_hvg, length(genes)))]
  sel <- sort(union(hvg, match(targets, genes)))
  cat("Genes kept:", length(sel), "(", n_hvg, "HVGs + target genes )\n")

  # Pass 2: read the selected genes block by block, stack at the end
  in_sel <- logical(length(genes))
  in_sel[sel] <- TRUE
  mats <- list()
  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    b_sel <- b[in_sel[b]]
    if (length(b_sel) == 0) next
    r <- read_block(min(b), max(b))
    if (is.null(r)) {
      mats[[length(mats) + 1]] <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0), x = numeric(0),
        dims = c(length(b_sel), n_keep))
    } else {
      local <- match(r$gene, b_sel)
      ok <- !is.na(local)
      mats[[length(mats) + 1]] <- Matrix::sparseMatrix(
        i = local[ok], j = r$cell[ok], x = r$value[ok],
        dims = c(length(b_sel), n_keep))
    }
    cat(sprintf("\rPass 2 (reading genes): block %d/%d", k, length(blocks)))
    rm(r); invisible(gc(verbose = FALSE))
  }
  cat("\nStacking blocks ...\n")
  expr <- do.call(rbind, mats)
  rm(mats); invisible(gc(verbose = FALSE))
  stopifnot(nrow(expr) == length(sel))
  dimnames(expr) <- list(genes[sel], rownames(obs)[keep])
  cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "cells,",
      format(length(expr@x), big.mark = ","), "non-zero values\n")

  list(expr = expr, perturbation = pert[keep], condition = cond[keep])
}

if (file.exists(cache)) {
  cat("Loading cached data:", cache, "\n")
  fr <- readRDS(cache)
} else {
  fr <- build_frangieh()
  cat("Saving cache ...\n")
  saveRDS(fr, cache)
}

# One pe_data object per condition, against that condition's own controls
make_condition <- function(cn) {
  cells <- fr$condition == cn
  pert <- fr$perturbation[cells]
  n <- table(pert)
  ok <- names(n)[n >= min_cells | names(n) == "control"]
  dropped <- setdiff(names(n), ok)
  if (length(dropped) > 0) {
    cat(cn, ": dropping", length(dropped), "perturbation(s) with fewer than",
        min_cells, "cells\n")
  }
  keep <- which(cells)[pert %in% ok]
  pe_data(fr$expr[, keep], fr$perturbation[keep], control = "control")
}
conditions <- c("monoculture", "IFNg", "coculture")
xs <- lapply(stats::setNames(conditions, conditions), make_condition)
for (cn in conditions) {
  cat("\n==", cn, "==\n")
  print(xs[[cn]])
}

# 2. Signal diagnostics per condition --------------------------------------------
sig_per <- list()
sig_sum <- list()
cat("\n== Signal diagnostics ==\n")
for (cn in conditions) {
  sig <- pe_signal(xs[[cn]], seed = 1)
  cat("\n--", cn, "\n")
  print(sig)
  sig_per[[cn]] <- cbind(condition = cn, sig$per_perturbation)
  sig_sum[[cn]] <- cbind(condition = cn, as.data.frame(sig$summary))
}
sig_per <- do.call(rbind, sig_per)
sig_sum <- do.call(rbind, sig_sum)
cat("\nMedians per condition:\n")
print(aggregate(cbind(effect_size, cosine_to_average, effective_genes,
                      reliability) ~ condition, data = sig_per,
                FUN = stats::median))
write.csv(sig_per, file.path(out_dir, "frangieh_signal_per_perturbation.csv"),
          row.names = FALSE)
write.csv(sig_sum, file.path(out_dir, "frangieh_signal_summary.csv"),
          row.names = FALSE)

# 3. Evaluate baselines over several splits, per condition -------------------------
raw <- list()
cal_raw <- list()
for (cn in conditions) {
  x <- xs[[cn]]
  for (sd in seeds) {
    s <- pe_split(x, type = "unseen_single", test_frac = 0.25, seed = sd)
    r <- pe_evaluate(list(baseline_control(x, s), baseline_mean(x, s),
                          baseline_linear(x, s)), x)
    cc <- pe_calibrate(r, pe_bounds(x, s, seed = sd))
    r$task <- cc$task <- cn
    r$seed <- cc$seed <- sd
    raw[[length(raw) + 1]] <- r
    cal_raw[[length(cal_raw) + 1]] <- cc
    cat(sprintf("%s: split %d/%d done\n", cn, sd, length(seeds)))
  }
}
res <- do.call(rbind, raw)
cal <- do.call(rbind, cal_raw)
write.csv(res, file.path(out_dir, "frangieh_per_perturbation.csv"),
          row.names = FALSE)
write.csv(cal, file.path(out_dir, "frangieh_calibrated_per_split.csv"),
          row.names = FALSE)

# 4. Summaries and paired comparisons ----------------------------------------------
keys <- unique(res[c("task", "method", "metric")])
summ <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  v <- res$value[res$task == keys$task[i] & res$method == keys$method[i] &
                   res$metric == keys$metric[i]]
  v <- v[!is.na(v)]
  data.frame(keys[i, ], mean = if (length(v)) mean(v) else NA_real_,
             se = if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_,
             n = length(v))
}))
rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "frangieh_summary.csv"), row.names = FALSE)

cmp <- do.call(rbind, lapply(unique(res$task), function(tk) {
  d <- res[res$task == tk, ]
  d$perturbation <- paste0("split", d$seed, ":", d$perturbation)
  out <- pe_compare(d[c("method", "perturbation", "metric", "value")],
                    reference = "mean", seed = 1)
  cbind(task = tk, out)
}))
write.csv(cmp, file.path(out_dir, "frangieh_compare_vs_mean.csv"),
          row.names = FALSE)

ck <- unique(cal[c("task", "method", "metric")])
cal_summ <- do.call(rbind, lapply(seq_len(nrow(ck)), function(i) {
  sel <- cal$task == ck$task[i] & cal$method == ck$method[i] &
    cal$metric == ck$metric[i]
  v <- cal$calibrated[sel]
  v <- v[!is.na(v)]
  data.frame(ck[i, ],
             lower = mean(cal$lower[sel], na.rm = TRUE),
             upper = mean(cal$upper[sel], na.rm = TRUE),
             calibrated = if (length(v)) mean(v) else NA_real_,
             sd_splits = if (length(v) > 1) stats::sd(v) else NA_real_)
}))
rownames(cal_summ) <- NULL
write.csv(cal_summ, file.path(out_dir, "frangieh_calibrated_summary.csv"),
          row.names = FALSE)

cat("\n== Mean over held-out perturbations (5 splits pooled) ==\n")
show <- summ
show$mean <- round(show$mean, 3)
show$se <- round(show$se, 3)
print(show[order(show$task, show$metric, show$method), ], row.names = FALSE)

cat("\n== Linear vs training mean (diff = linear - mean) ==\n")
show_cmp <- cmp[cmp$method_a == "linear",
                c("task", "metric", "n", "diff", "ci_low", "ci_high",
                  "p_adj", "favours")]
show_cmp[, c("diff", "ci_low", "ci_high")] <-
  round(show_cmp[, c("diff", "ci_low", "ci_high")], 3)
show_cmp$p_adj <- signif(show_cmp$p_adj, 2)
print(show_cmp[order(show_cmp$task, show_cmp$metric), ], row.names = FALSE)

cat("\n== Calibrated scores of the linear model (0 = best uninformative baseline, 1 = replicate) ==\n")
show_cal <- cal_summ[cal_summ$method == "linear", ]
num <- c("lower", "upper", "calibrated", "sd_splits")
show_cal[num] <- round(show_cal[num], 3)
print(show_cal[order(show_cal$task, show_cal$metric), ], row.names = FALSE)

cat("\nResults written to", out_dir, "\n")
