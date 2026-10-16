# Figure 3 (draft): calibrated scores on Norman et al. 2019
# ---------------------------------------------------------------------------
# Reads analysis/output/norman_calibrated_per_split.csv written by
# real02_norman_baselines.R (run that first) and plots, for each task and
# metric, every method's calibrated score:
#   0 = the better uninformative baseline (no change / training mean)
#   1 = split-half technical replicate, corrected to full sample size
# Points are means over the 5 random splits; bars are +/- 1 SD across splits.
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real03_calibration_figure.R"))

# Package folder on this computer's Desktop (works for any Windows user name).
# If your Desktop is elsewhere (e.g. OneDrive), set pkg to the full path.
pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
out_dir <- file.path(pkg, "analysis", "output")
in_file <- file.path(out_dir, "norman_calibrated_per_split.csv")
if (!file.exists(in_file)) {
  stop("Run real02_norman_baselines.R first: ", in_file, " not found.")
}
cal <- read.csv(in_file, stringsAsFactors = FALSE)

# 1. Mean and SD over splits ------------------------------------------------------
keys <- unique(cal[c("task", "method", "metric")])
fig <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  v <- cal$calibrated[cal$task == keys$task[i] &
                        cal$method == keys$method[i] &
                        cal$metric == keys$metric[i]]
  v <- v[!is.na(v)]
  data.frame(keys[i, ],
             value = if (length(v)) mean(v) else NA_real_,
             sd = if (length(v) > 1) stats::sd(v) else NA_real_,
             n_splits = length(v))
}))
rownames(fig) <- NULL
write.csv(fig, file.path(out_dir, "fig3_calibrated_summary.csv"),
          row.names = FALSE)

# 2. Plot -------------------------------------------------------------------------
tasks <- c(unseen_single = "A  Unseen single perturbations",
           unseen_combo = "B  Unseen pairs")
metrics <- c(pds = "Discrimination (PDS)",
             pearson_delta_pert = "PearsonDelta vs centroid",
             pearson_delta_top = "PearsonDelta, top 20 genes",
             pearson_delta = "PearsonDelta vs control",
             mse_top = "MSE, top 20 genes",
             mse = "MSE")
labels <- c(control = "No change", mean = "Training mean",
            linear = "Linear", additive = "Additive")
cols <- c(control = "#2a78d6", mean = "#eb6834", linear = "#1baf7a",
          additive = "#eda100")
pchs <- c(control = 16, mean = 17, linear = 15, additive = 18)
offset <- c(control = -0.22, mean = 0, linear = 0.22, additive = 0.22)
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"
surface <- "#fcfcfb"

xr <- range(c(fig$value - fig$sd, fig$value + fig$sd, 0, 1), na.rm = TRUE)
xr <- c(min(-0.1, xr[1]), max(1.1, xr[2]))
xr <- xr + c(-1, 1) * 0.03 * diff(xr)

draw_fig3 <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(oma = c(1.6, 0, 0, 0))
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE),
         heights = c(0.12, 0.88))
  par(bg = surface, col.axis = ink2, col.lab = ink2, fg = ink2)

  # Legend row
  par(mar = c(0, 0, 0, 0))
  plot.new()
  # Equal slot width for every entry so the spacing is even
  slot <- max(strwidth(labels, cex = 1.05)) * 1.15
  legend("center", legend = labels, col = cols[names(labels)],
         pch = pchs[names(labels)], lwd = 2, horiz = TRUE, bty = "n",
         text.col = ink, cex = 1.05, x.intersp = 0.6, seg.len = 1.4,
         text.width = slot)

  for (tk in names(tasks)) {
    d <- fig[fig$task == tk, ]
    ny <- length(metrics)
    par(mar = c(4.6, if (tk == "unseen_single") 13 else 2, 3, 1.2))
    plot(NA, xlim = xr, ylim = c(0.4, ny + 0.6), axes = FALSE,
         xlab = "Calibrated score", ylab = "")
    xt <- pretty(xr)
    abline(v = xt, col = grid_col)
    abline(h = seq_len(ny - 1) + 0.5, col = grid_col, lty = 1)
    abline(v = c(0, 1), col = ink2, lty = 2)
    axis(1, at = xt, lwd = 0, lwd.ticks = 1, col.ticks = grid_col)
    if (tk == "unseen_single") {
      axis(2, at = seq_len(ny), labels = metrics, las = 1, lwd = 0,
           cex.axis = 0.95, col.axis = ink)
    }
    mtext(c("uninformative", "replicate"), side = 3, at = c(0, 1),
          line = 0.1, cex = 0.75, col = ink2)

    for (k in seq_len(ny)) {
      mt <- names(metrics)[k]
      for (m in intersect(names(labels), d$method)) {
        row <- d[d$metric == mt & d$method == m, ]
        if (nrow(row) == 0 || is.na(row$value)) next
        y <- k + offset[[m]]
        if (!is.na(row$sd) && row$sd > 0) {
          arrows(row$value - row$sd, y, row$value + row$sd, y,
                 angle = 90, code = 3, length = 0.03, col = cols[[m]],
                 lwd = 1.5)
        }
        points(row$value, y, pch = pchs[[m]], col = cols[[m]], cex = 1.5)
      }
    }
    title(main = tasks[[tk]], adj = 0, col.main = ink, cex.main = 1.1,
          font.main = 2, line = 1.4)
  }
  mtext(paste("Points: mean over 5 random splits; bars: +/- 1 SD across",
              "splits. A missing point means the metric is undefined for",
              "that baseline (it predicts no change from the metric's",
              "reference)."),
        side = 1, outer = TRUE, line = 0.5, adj = 0, at = 0.02, cex = 0.8,
        col = ink2)
}

png(file.path(out_dir, "fig3_calibrated_scores.png"), width = 3600,
    height = 1900, res = 300)
draw_fig3()
dev.off()
pdf(file.path(out_dir, "fig3_calibrated_scores.pdf"), width = 12,
    height = 6.3)
draw_fig3()
dev.off()

show <- fig
show$value <- round(show$value, 3)
show$sd <- round(show$sd, 3)
print(show[order(show$task, show$metric, show$method),
           c("task", "metric", "method", "value", "sd")], row.names = FALSE)
message("Figure 3 written to ", out_dir)
