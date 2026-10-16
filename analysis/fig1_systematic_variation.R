# Figure 1 (draft): systematic variation flatters uninformative baselines
# ---------------------------------------------------------------------------
# Question: when every perturbed cell shares a common shift ("systematic
# variation"), do common metrics still reward models for knowing the
# perturbation-specific effects?
#
# Design: simulate datasets with increasing shared shift (`shared_effect`),
# hold out unseen single perturbations, and score four predictions:
#   * No change         - baseline_control(): knows nothing
#   * Training mean     - baseline_mean(): knows only the shared shift
#   * Linear            - baseline_linear(): Ahlmann-Eltze et al. 2025
#   * True effects only - control + the exact simulated perturbation-specific
#                         effects, but no shared shift: knows everything
#                         perturbation-specific, nothing systematic
# Metrics: PearsonDelta vs control, MSE, perturbation discrimination score.
#
# Expected (checked with a quick numpy prototype of the same simulation):
#   * PearsonDelta and MSE flip from about shared_effect = 0.5 on: the training
#     mean, which knows nothing perturbation-specific, beats the model that
#     knows the true specific effects (prototype at 0.5: r 0.81 vs 0.55,
#     MSE 0.13 vs 0.26).
#   * The discrimination score keeps the training mean at chance (0.5) and the
#     true-effects model near 1; with the default L1 distance it starts to
#     slip only at very strong shifts (about 0.89 at 3 SD in the prototype).
#
# Section 4 adds a paired statistical comparison (pe_compare) of the
# true-effects model against the training mean at every shift level.
#
# Run from R (no installation needed; uses the package source directly):
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "fig1_systematic_variation.R"))
# Outputs go to analysis/output/.

# Package folder on this computer's Desktop (works for any Windows user name).
# If your Desktop is elsewhere (e.g. OneDrive), set pkg to the full path.
pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
devtools::load_all(pkg, quiet = TRUE)

out_dir <- file.path(pkg, "analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

shared_levels <- c(0, 0.25, 0.5, 1, 1.5, 2, 3)
replicates <- 1:5
metrics <- c("pearson_delta", "mse", "pds")

# A prediction that knows the true perturbation-specific effects but not the
# shift shared by all perturbed cells.
true_effects_only <- function(x, s) {
  ctrl <- pe_pseudobulk(x, x$control)[, 1]
  pred <- ctrl + x$sim$effects[names(ctrl), s$test, drop = FALSE]
  pe_prediction(pred, s, method = "true_effects_only")
}

# 1. Simulate and evaluate ------------------------------------------------------

rows <- list()
raw <- list()   # per-perturbation results, for pe_compare() below
for (lev in shared_levels) {
  for (r in replicates) {
    x <- pe_simulate(shared_effect = lev, effect_sd = 1, seed = r)
    s <- pe_split(x, type = "unseen_single", test_frac = 0.3, seed = r)
    preds <- list(
      baseline_control(x, s),
      baseline_mean(x, s),
      baseline_linear(x, s, targets = x$sim$targets),
      true_effects_only(x, s)
    )
    res <- pe_evaluate(preds, x, metrics = metrics)
    sm <- pe_summary(res)
    sm$shared_effect <- lev
    sm$replicate <- r
    rows[[length(rows) + 1]] <- sm
    # Label perturbations by replicate so replicates can be pooled.
    res$perturbation <- paste0("rep", r, ":", res$perturbation)
    res$shared_effect <- lev
    raw[[length(raw) + 1]] <- res
  }
}
per_rep <- do.call(rbind, rows)
write.csv(per_rep, file.path(out_dir, "fig1_per_replicate.csv"),
          row.names = FALSE)

# 2. Average over replicates ---------------------------------------------------

keys <- unique(per_rep[c("method", "metric", "shared_effect")])
agg <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  v <- per_rep$mean[per_rep$method == keys$method[i] &
                      per_rep$metric == keys$metric[i] &
                      per_rep$shared_effect == keys$shared_effect[i]]
  v <- v[!is.na(v)]
  data.frame(keys[i, ],
             value = if (length(v)) mean(v) else NA_real_,
             sd = if (length(v) > 1) stats::sd(v) else NA_real_,
             n_rep = length(v))
}))
rownames(agg) <- NULL
write.csv(agg, file.path(out_dir, "fig1_summary.csv"), row.names = FALSE)

# 3. Plot --------------------------------------------------------------------
# Fixed categorical order (validated: CVD-safe; two hues are light on white,
# so every series also has its own point shape and the data ship as CSV).

methods <- c("control", "mean", "linear", "true_effects_only")
labels <- c(control = "No change (control)",
            mean = "Training mean",
            linear = "Linear (Ahlmann-Eltze)",
            true_effects_only = "True specific effects only")
cols <- c(control = "#2a78d6", mean = "#eb6834",
          linear = "#1baf7a", true_effects_only = "#eda100")
pchs <- c(control = 16, mean = 17, linear = 15, true_effects_only = 18)
# Mean and linear baselines nearly coincide in this simulation (its effects
# carry no co-expression structure for the linear model to learn), as do the
# no-change and true-effects MSE at large shifts. A dashed line for the mean
# and a small horizontal dodge keep every series visible.
ltys <- c(control = 1, mean = 2, linear = 1, true_effects_only = 1)
dodge <- c(control = -0.045, mean = -0.015, linear = 0.015,
           true_effects_only = 0.045)
tick_label <- function(v) format(v, scientific = FALSE, drop0trailing = TRUE,
                                 trim = TRUE)

panels <- list(
  list(metric = "pearson_delta", title = "A  PearsonDelta vs control",
       ylab = "Pearson r (higher is better)", ref = NA, log = FALSE),
  list(metric = "mse", title = "B  Mean squared error",
       ylab = "MSE, log scale (lower is better)", ref = NA, log = TRUE),
  list(metric = "pds", title = "C  Discrimination score (L1)",
       ylab = "PDS (higher is better)", ref = 0.5, log = FALSE)
)

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; surface <- "#fcfcfb"

draw_figure <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE),
         heights = c(0.13, 0.87))
  par(bg = surface, col.axis = ink2, col.lab = ink2, fg = ink2,
      family = "sans")

  # Shared legend row
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = labels[methods], col = cols[methods],
         pch = pchs[methods], lty = ltys[methods], lwd = 2, horiz = TRUE,
         bty = "n", cex = 1.05, text.col = ink, x.intersp = 0.6,
         seg.len = 2.2)

  for (pn in panels) {
    d <- agg[agg$metric == pn$metric, ]
    if (pn$log) {
      lim <- range(d$value, na.rm = TRUE) * c(0.6, 1.6)
      d$lo <- pmax(d$value - d$sd, lim[1])
    } else {
      lim <- range(c(d$value - d$sd, d$value + d$sd, pn$ref), na.rm = TRUE)
      lim <- lim + c(-1, 1) * 0.04 * diff(lim)
      d$lo <- d$value - d$sd
    }
    d$hi <- d$value + d$sd
    par(mar = c(5, 5.2, 2.6, 1), mgp = c(3.4, 0.7, 0))
    plot(NA, xlim = range(shared_levels) + c(-0.08, 0.08), ylim = lim,
         axes = FALSE, log = if (pn$log) "y" else "",
         xlab = "Shared shift across perturbed cells (SD)", ylab = pn$ylab)
    yt <- axTicks(2)
    abline(h = yt, col = grid_col, lwd = 1)
    axis(1, at = shared_levels, labels = tick_label(shared_levels),
         lwd = 0, lwd.ticks = 1, col.ticks = grid_col, gap.axis = 0,
         mgp = c(3, 0.7, 0))
    axis(2, at = yt, labels = tick_label(yt), las = 1, lwd = 0,
         lwd.ticks = 1, col.ticks = grid_col)
    if (!is.na(pn$ref)) {
      abline(h = pn$ref, lty = 3, col = ink2)
      text(max(shared_levels) + 0.08, pn$ref, "chance", adj = c(1, 1.6),
           cex = 0.85, col = ink2)
    }
    for (m in methods) {
      dm <- d[d$method == m & !is.na(d$value), ]
      if (nrow(dm) == 0) next   # e.g. no-change baseline has no PearsonDelta
      dm <- dm[order(dm$shared_effect), ]
      xs <- dm$shared_effect + dodge[m]
      ok <- !is.na(dm$sd) & dm$hi > dm$lo
      if (any(ok)) {
        arrows(xs[ok], dm$lo[ok], xs[ok], dm$hi[ok],
               angle = 90, code = 3, length = 0.03, col = cols[m], lwd = 1)
      }
      lines(xs, dm$value, col = cols[m], lwd = 2, lty = ltys[m])
      points(xs, dm$value, col = cols[m], pch = pchs[m], cex = 1.3)
    }
    title(main = pn$title, adj = 0, font.main = 2, col.main = ink,
          cex.main = 1.1)
  }
}

png(file.path(out_dir, "fig1_systematic_variation.png"),
    width = 3600, height = 1350, res = 300)
draw_figure()
dev.off()

pdf(file.path(out_dir, "fig1_systematic_variation.pdf"),
    width = 12, height = 4.5)
draw_figure()
dev.off()

# 4. Paired comparison: true specific effects vs training mean ------------------
# Pools the held-out perturbations of all replicates at each shift level
# (5 x 6 = 30 paired observations) and asks, per metric, which of the two
# predictions is better, with a bootstrap 95% CI for the mean difference.

per_pert <- do.call(rbind, raw)
cmp <- do.call(rbind, lapply(shared_levels, function(lev) {
  d <- per_pert[per_pert$shared_effect == lev &
                  per_pert$method %in% c("true_effects_only", "mean"), ]
  d$shared_effect <- NULL
  out <- pe_compare(d, reference = "mean", n_boot = 2000, seed = 1)
  cbind(shared_effect = lev, out)
}))
write.csv(cmp, file.path(out_dir, "fig1_compare_true_vs_mean.csv"),
          row.names = FALSE)

cat("\nTrue specific effects (a) vs training mean (b), pooled replicates:\n")
show_cmp <- cmp[, c("shared_effect", "metric", "diff", "ci_low", "ci_high",
                    "p_adj", "favours")]
show_cmp[, c("diff", "ci_low", "ci_high")] <-
  round(show_cmp[, c("diff", "ci_low", "ci_high")], 3)
show_cmp$p_adj <- signif(show_cmp$p_adj, 2)
print(show_cmp[order(show_cmp$metric, show_cmp$shared_effect), ],
      row.names = FALSE)

# 5. Console summary at the extremes ------------------------------------------

show <- agg[agg$shared_effect %in% range(shared_levels),
            c("shared_effect", "metric", "method", "value")]
show$value <- round(show$value, 3)
show <- show[order(show$metric, show$shared_effect, show$method), ]
print(show, row.names = FALSE)
message("Figure and tables written to ", out_dir)
