# Real data, step 2: baselines and metrics on Norman et al. 2019
# ---------------------------------------------------------------------------
# Data: scPerturb NormanWeissman2019_filtered.h5ad (Zenodo 10044268, CC BY 4.0)
#   111,445 K562 cells; obs$perturbation: "control", single genes ("KLF1"),
#   pairs joined by "_" ("CEBPE_RUNX1T1"); obs$ncounts = total UMI per cell;
#   X is stored by gene (CSC), raw counts.
#
# Memory: the full matrix would need ~4 GB in R, so it is never loaded whole.
#   Pass 1 reads X in gene blocks to find highly variable genes; pass 2 reads
#   only the selected genes. Cells are capped per perturbation. The processed
#   dataset is cached as analysis/data/norman_pe_data.rds; later runs start
#   from the cache (delete it to rebuild with different settings).
#
# Analysis: 5 random splits for each task
#   * unseen single perturbations: no change, training mean, linear model
#   * unseen combinations:         no change, training mean, additive model
# Outputs (analysis/output/): per-perturbation results, summary, paired
# comparisons against the training mean, calibrated scores (pe_bounds +
# pe_calibrate), signal diagnostics (pe_signal), and a figure.
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real02_norman_baselines.R"))

# Package folder on this computer's Desktop (works for any Windows user name).
# If your Desktop is elsewhere (e.g. OneDrive), set pkg to the full path.
pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
devtools::load_all(pkg, quiet = TRUE)

data_dir <- file.path(pkg, "analysis", "data")
out_dir <- file.path(pkg, "analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
h5ad <- file.path(data_dir, "NormanWeissman2019_filtered.h5ad")
cache <- file.path(data_dir, "norman_pe_data.rds")

# Settings
n_hvg <- 5000          # highly variable genes kept (plus all target genes)
max_cells_pert <- 300  # cells kept per perturbation
max_cells_ctrl <- 3000 # control cells kept
block_genes <- 1500    # genes read from disk at a time
seeds <- 1:5           # random splits per task

# 1. Build (or load) the pe_data object ----------------------------------------

build_norman <- function() {
  enc <- rhdf5::h5readAttributes(h5ad, "/X")[["encoding-type"]]
  cat("X encoding:", enc, "\n")
  if (!identical(enc, "csc_matrix")) {
    stop("Expected X stored as csc_matrix (by gene); found ", enc, ".")
  }

  ad <- anndataR::read_h5ad(h5ad, as = "HDF5AnnData")
  obs <- ad$obs
  genes <- make.unique(rownames(ad$var))
  pert <- as.character(obs$perturbation)
  libsize <- obs$ncounts
  n_cells <- length(pert)
  # Close anndataR's file handle so the block reads below do not warn about
  # an already-open HDF5 file.
  rm(ad)
  rhdf5::h5closeAll()

  # Label sanity checks against nperts
  stopifnot(all(pert[obs$nperts == 0] == "control"),
            !any(grepl("_", pert[obs$nperts == 1], fixed = TRUE)),
            all(grepl("_", pert[obs$nperts == 2], fixed = TRUE)))

  # Drop perturbations whose target gene is not among the measured genes
  labels <- setdiff(unique(pert), "control")
  parts <- strsplit(labels, "_", fixed = TRUE)
  absent <- labels[vapply(parts, function(g) !all(g %in% genes), logical(1))]
  if (length(absent) > 0) {
    cat("Dropping", length(absent), "perturbation(s) with unmeasured targets:",
        paste(absent, collapse = ", "), "\n")
  }
  targets <- unique(unlist(parts))
  targets <- targets[targets %in% genes]

  # Cap cells per perturbation
  set.seed(1)
  by_label <- split(seq_len(n_cells), pert)
  by_label <- by_label[!names(by_label) %in% absent]
  keep <- sort(unlist(lapply(names(by_label), function(l) {
    idx <- by_label[[l]]
    cap <- if (l == "control") max_cells_ctrl else max_cells_pert
    if (length(idx) > cap) idx[sample.int(length(idx), cap)] else idx
  }), use.names = FALSE))
  n_keep <- length(keep)
  cell_map <- integer(n_cells)
  cell_map[keep] <- seq_len(n_keep)
  scale <- 1e4 / libsize
  cat("Cells kept:", n_keep, "of", n_cells, "\n")

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

  # Pass 2: read the selected genes. Each block becomes a small sparse matrix
  # holding only its own selected genes; the blocks are stacked at the end.
  in_sel <- logical(length(genes))
  in_sel[sel] <- TRUE
  mats <- list()
  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    b_sel <- b[in_sel[b]]          # selected genes in this block, in order
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

  pe_data(expr, pert[keep], control = "control", combo_sep = "_")
}

if (file.exists(cache)) {
  cat("Loading cached dataset:", cache, "\n")
  x <- readRDS(cache)
} else {
  x <- build_norman()
  cat("Saving cache ...\n")
  saveRDS(x, cache)
}
print(x)

# 2. Signal diagnostics -------------------------------------------------------------
# pe_signal(): cosine similarity between each perturbation's change from control
# and the average change over all perturbations (the Systema statistic; near 1 =
# perturbations mostly share one common shift), plus effect size, effective
# number of changed genes and split-half reliability.
sig <- pe_signal(x, seed = 1)
print(sig)
cat("\nBy perturbation type (medians):\n")
sp <- sig$per_perturbation
sp$type <- ifelse(grepl(x$combo_sep, sp$perturbation, fixed = TRUE),
                  "combo", "single")
print(aggregate(cbind(effect_size, cosine_to_average, effective_genes,
                      reliability) ~ type, data = sp, FUN = stats::median))
write.csv(sp, file.path(out_dir, "norman_systematic_variation.csv"),
          row.names = FALSE)
write.csv(as.data.frame(sig$summary),
          file.path(out_dir, "norman_signal_summary.csv"), row.names = FALSE)

# 3. Evaluate baselines over several splits ---------------------------------------
raw <- list()
cal_raw <- list()   # calibrated scores per split (pe_bounds + pe_calibrate)
for (sd in seeds) {
  s1 <- pe_split(x, type = "unseen_single", test_frac = 0.25, seed = sd)
  r1 <- pe_evaluate(list(baseline_control(x, s1), baseline_mean(x, s1),
                         baseline_linear(x, s1)), x)
  c1 <- pe_calibrate(r1, pe_bounds(x, s1, seed = sd))
  r1$task <- c1$task <- "unseen_single"

  s2 <- pe_split(x, type = "unseen_combo", test_frac = 0.25, seed = sd)
  r2 <- pe_evaluate(list(baseline_control(x, s2), baseline_mean(x, s2),
                         baseline_additive(x, s2)), x)
  c2 <- pe_calibrate(r2, pe_bounds(x, s2, seed = sd))
  r2$task <- c2$task <- "unseen_combo"

  r <- rbind(r1, r2)
  r$seed <- sd
  raw[[length(raw) + 1]] <- r
  cc <- rbind(c1, c2)
  cc$seed <- sd
  cal_raw[[length(cal_raw) + 1]] <- cc
  cat(sprintf("Split %d/%d done\n", sd, length(seeds)))
}
res <- do.call(rbind, raw)
cal <- do.call(rbind, cal_raw)
write.csv(cal, file.path(out_dir, "norman_calibrated_per_split.csv"),
          row.names = FALSE)
write.csv(res, file.path(out_dir, "norman_per_perturbation.csv"),
          row.names = FALSE)

# 4. Summaries and paired comparisons ---------------------------------------------
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
write.csv(summ, file.path(out_dir, "norman_summary.csv"), row.names = FALSE)

cmp <- do.call(rbind, lapply(unique(res$task), function(tk) {
  d <- res[res$task == tk, ]
  # Pool splits: label perturbations by split so each row is unique
  d$perturbation <- paste0("split", d$seed, ":", d$perturbation)
  out <- pe_compare(d[c("method", "perturbation", "metric", "value")],
                    reference = "mean", seed = 1)
  cbind(task = tk, out)
}))
write.csv(cmp, file.path(out_dir, "norman_compare_vs_mean.csv"),
          row.names = FALSE)

cat("\n== Mean over held-out perturbations (5 splits pooled) ==\n")
show <- summ
show$mean <- round(show$mean, 3)
show$se <- round(show$se, 3)
print(show[order(show$task, show$metric, show$method), ], row.names = FALSE)

cat("\n== Each method vs training mean (diff = method - mean) ==\n")
show_cmp <- cmp[, c("task", "metric", "method_a", "n", "diff", "ci_low",
                    "ci_high", "favours")]
show_cmp[, c("diff", "ci_low", "ci_high")] <-
  round(show_cmp[, c("diff", "ci_low", "ci_high")], 3)
print(show_cmp[order(show_cmp$task, show_cmp$metric), ], row.names = FALSE)

# Calibrated scores: 0 = best uninformative baseline, 1 = split-half replicate.
# Averaged over the 5 splits (sd across splits in brackets).
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
write.csv(cal_summ, file.path(out_dir, "norman_calibrated_summary.csv"),
          row.names = FALSE)

cat("\n== Calibrated scores (0 = best uninformative baseline, 1 = replicate) ==\n")
show_cal <- cal_summ
num <- c("lower", "upper", "calibrated", "sd_splits")
show_cal[num] <- round(show_cal[num], 3)
print(show_cal[order(show_cal$task, show_cal$metric, show_cal$method), ],
      row.names = FALSE)

# 5. Figure -----------------------------------------------------------------------
# Same categorical order as Figure 1; additive takes slot 4.
labels <- c(control = "No\nchange", mean = "Training\nmean",
            linear = "Linear", additive = "Additive")
tick_label <- function(v) format(v, scientific = FALSE, drop0trailing = TRUE,
                                 trim = TRUE)
cols <- c(control = "#2a78d6", mean = "#eb6834", linear = "#1baf7a",
          additive = "#eda100")
pchs <- c(control = 16, mean = 17, linear = 15, additive = 18)
tasks <- c(unseen_single = "Unseen singles", unseen_combo = "Unseen pairs")
mets <- c(pearson_delta = "PearsonDelta vs control",
          pearson_delta_pert = "PearsonDelta vs centroid",
          mse = "MSE (lower is better)",
          pds = "Discrimination score")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; surface <- "#fcfcfb"

draw_real <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(2, 4), mar = c(3.6, 4.4, 3, 1), oma = c(2, 0, 0, 0),
      bg = surface, col.axis = ink2, col.lab = ink2, fg = ink2)
  for (tk in names(tasks)) {
    for (mt in names(mets)) {
      d <- summ[summ$task == tk & summ$metric == mt, ]
      ms <- intersect(names(labels), d$method)
      d <- d[match(ms, d$method), ]
      lo <- d$mean - 1.96 * d$se
      hi <- d$mean + 1.96 * d$se
      rng <- range(c(lo, hi, d$mean, if (mt == "pds") 0.5), na.rm = TRUE)
      if (!all(is.finite(rng))) rng <- c(0, 1)
      rng <- rng + c(-1, 1) * 0.08 * max(diff(rng), 1e-6)
      plot(NA, xlim = c(0.5, length(ms) + 0.5), ylim = rng, axes = FALSE,
           xlab = "", ylab = "")
      yt <- axTicks(2)
      abline(h = yt, col = grid_col)
      axis(2, at = yt, labels = tick_label(yt), las = 1, lwd = 0,
           lwd.ticks = 1, col.ticks = grid_col)
      axis(1, at = seq_along(ms), labels = labels[ms], las = 1, lwd = 0,
           padj = 0.6, cex.axis = 0.95)
      if (mt == "pds") abline(h = 0.5, lty = 3, col = ink2)
      for (j in seq_along(ms)) {
        m <- ms[j]
        if (is.na(d$mean[j])) {
          text(j, mean(rng), "n/a", col = ink2, cex = 0.9)
          next
        }
        if (!is.na(d$se[j]) && hi[j] > lo[j]) {
          arrows(j, lo[j], j, hi[j], angle = 90, code = 3, length = 0.04,
                 col = cols[m], lwd = 1.5)
        }
        points(j, d$mean[j], pch = pchs[m], col = cols[m], cex = 1.6)
      }
      title(main = paste0(tasks[tk], ": ", mets[mt]), adj = 0,
            col.main = ink, cex.main = 1, font.main = 2)
    }
  }
  mtext(paste("Points: mean over held-out perturbations (5 random splits",
              "pooled); bars: 95% CI. n/a: metric undefined for that",
              "baseline. Paired differences vs the training mean are in",
              "norman_compare_vs_mean.csv."),
        side = 1, outer = TRUE, line = 0.6, adj = 0, cex = 0.8, col = ink2)
}

png(file.path(out_dir, "norman_baselines.png"), width = 4200, height = 2400,
    res = 300)
draw_real()
dev.off()
pdf(file.path(out_dir, "norman_baselines.pdf"), width = 14, height = 8)
draw_real()
dev.off()

message("Results written to ", out_dir)
