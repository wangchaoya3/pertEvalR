# Figure 4: dataset signal, headroom and reliability-stratified comparison
# ---------------------------------------------------------------------------
# Reads the CSV files written by real02 (Norman), real05 and real06
# (Frangieh); no data are recomputed.
#   (a) per-perturbation reliability vs cosine to the average shift
#   (b) headroom: best uninformative baseline -> technical replicate, with the
#       model score, for PearsonDelta
#   (c) linear minus training mean (95% CI) by reliability stratum
#   (d) PDS of an uninformative prediction (training mean) vs effect size
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real07_fig4_signal.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
out_dir <- file.path(pkg, "analysis", "output")
rd <- function(f) {
  p <- file.path(out_dir, f)
  if (!file.exists(p)) stop("Missing ", f, "; run real02, real05 and real06 first.")
  utils::read.csv(p, stringsAsFactors = FALSE)
}

nor_sig <- rd("norman_systematic_variation.csv")
nor_res <- rd("norman_per_perturbation.csv")
nor_cal <- rd("norman_calibrated_summary.csv")
fr_sig <- rd("frangieh_signal_per_perturbation.csv")
fr_res <- rd("frangieh_reliable_per_perturbation.csv")
fr_cmp <- rd("frangieh_reliable_compare_vs_mean.csv")
fr_cal <- rd("frangieh_reliable_calibrated.csv")

if (!"reliability" %in% names(nor_sig)) {
  stop("norman_systematic_variation.csv has no reliability column; ",
       "re-run real02 with pertEvalR >= 0.0.0.9005.")
}

# Style (same tokens as Figures 1 and 3) --------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; surface <- "#fcfcfb"
ds_cols <- c(norman = "#2a78d6", monoculture = "#1baf7a", IFNg = "#eb6834",
             coculture = "#eda100")
ds_labels <- c(norman = "Norman (K562)", monoculture = "Frangieh monoculture",
               IFNg = "Frangieh IFN-\u03b3", coculture = "Frangieh co-culture")
ifn_genes <- c("IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1")
min_rel <- 0.5
tick_label <- function(v) format(v, scientific = FALSE, drop0trailing = TRUE,
                                 trim = TRUE)
alpha <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)
panel_label <- function(lab) {
  mtext(lab, side = 3, line = 1.2, adj = 0,
        at = graphics::grconvertX(0.015, "nfc", "user"), font = 2,
        cex = 1.25, col = ink)
}
cs <- 0.9   # small text (legends, annotations, axis labels in b and c)
grid_y <- function(at) abline(h = at, col = grid_col)
grid_x <- function(at) abline(v = at, col = grid_col)

# Per-perturbation signal table for all four datasets
sig <- rbind(
  data.frame(dataset = "norman", perturbation = nor_sig$perturbation,
             cosine = nor_sig$cosine_to_average,
             reliability = nor_sig$reliability,
             effect_size = nor_sig$effect_size),
  data.frame(dataset = fr_sig$condition, perturbation = fr_sig$perturbation,
             cosine = fr_sig$cosine_to_average,
             reliability = fr_sig$reliability,
             effect_size = fr_sig$effect_size)
)

# (a) Reliability vs systematic variation ---------------------------------------
panel_a <- function() {
  par(mar = c(3.6, 3.8, 2.6, 0.8))
  x_min <- min(-0.25, floor(min(sig$cosine, na.rm = TRUE) * 4) / 4)
  xt <- seq(x_min, 1, 0.25)
  plot(NA, xlim = c(x_min, 1), ylim = c(-0.2, 1), axes = FALSE,
       xlab = "Cosine to the average shift",
       ylab = "Split-half reliability")
  grid_x(xt); grid_y(seq(0, 1, 0.25))
  axis(1, at = xt, labels = tick_label(xt),
       lwd = 0, lwd.ticks = 1, col.ticks = grid_col)
  axis(2, at = seq(0, 1, 0.25), labels = tick_label(seq(0, 1, 0.25)),
       las = 1, lwd = 0, lwd.ticks = 1, col.ticks = grid_col)
  abline(h = min_rel, lty = 3, col = ink2)
  abline(v = 0, col = ink2, lwd = 0.8)
  text(x_min, min_rel, paste("reliability", min_rel), adj = c(0, -0.4),
       cex = cs, col = ink2)
  for (d in c("monoculture", "IFNg", "coculture", "norman")) {
    s <- sig[sig$dataset == d, ]
    points(s$cosine, pmax(s$reliability, -0.2), pch = 16,
           cex = 0.55, col = alpha(ds_cols[d], 0.55))
  }
  ifn <- sig[sig$dataset %in% c("IFNg", "coculture") &
               sig$perturbation %in% ifn_genes, ]
  points(ifn$cosine, ifn$reliability, pch = 21, cex = 0.8, lwd = 0.7,
         col = ink2, bg = ds_cols[ifn$dataset])
  # Legend in the empty lower-left corner; the IFN-gamma pathway knockouts
  # get their own entry instead of a text label over the point cloud.
  legend("bottomleft", legend = c(ds_labels, "IFN-\u03b3 pathway KO"),
         pch = c(rep(16, length(ds_labels)), 21),
         col = c(ds_cols, ink2), pt.bg = c(ds_cols, "#bdbcb8"),
         pt.cex = c(rep(1, length(ds_labels)), 0.8), pt.lwd = 0.7,
         bty = "n", cex = cs, text.col = ink2)
  title(main = "Signal per perturbation", adj = 0, col.main = ink,
        cex.main = 1, font.main = 2)
  panel_label("a")
}

# (b) Headroom for PearsonDelta ----------------------------------------------------
headroom_rows <- function() {
  rows <- list()
  add <- function(label, ds, lower, upper, score, model) {
    rows[[length(rows) + 1]] <<- data.frame(label = label, dataset = ds,
                                            lower = lower, upper = upper,
                                            score = score, model = model)
  }
  for (tk in c("unseen_combo", "unseen_single")) {
    m <- if (tk == "unseen_combo") "additive" else "linear"
    r <- nor_cal[nor_cal$task == tk & nor_cal$method == m &
                   nor_cal$metric == "pearson_delta", ]
    if (nrow(r) == 1) {
      add(if (tk == "unseen_combo") "Norman pairs" else "Norman singles",
          "norman", r$lower, r$upper,
          r$lower + r$calibrated * (r$upper - r$lower), m)
    }
  }
  for (cn in c("IFNg", "coculture", "monoculture")) {
    for (st in c("reliable", "unreliable")) {
      r <- fr_cal[fr_cal$condition == cn & fr_cal$stratum == st &
                    fr_cal$method == "linear" &
                    fr_cal$metric == "pearson_delta", ]
      if (nrow(r) == 1) {
        add(sprintf("%s\n%s (%d)", sub("Frangieh ", "", ds_labels[cn]), st,
                    r$n_distinct), cn, r$lower, r$upper, r$score, "linear")
      }
    }
  }
  do.call(rbind, rows)
}

panel_b <- function() {
  h <- headroom_rows()
  n <- nrow(h)
  y <- rev(seq_len(n))
  par(mar = c(3.6, 7.6, 2.6, 2.6))
  plot(NA, xlim = c(0, 1), ylim = c(0.5, n + 0.5), axes = FALSE,
       xlab = "Pearson\u0394 vs control", ylab = "")
  grid_x(seq(0, 1, 0.25))
  axis(1, at = seq(0, 1, 0.25), labels = tick_label(seq(0, 1, 0.25)),
       lwd = 0, lwd.ticks = 1, col.ticks = grid_col)
  axis(2, at = y, labels = h$label, las = 1, lwd = 0, cex.axis = cs,
       padj = 0.5)
  segments(h$lower, y, h$upper, y, lwd = 6, col = alpha(ds_cols[h$dataset], 0.25),
           lend = 1)
  segments(h$lower, y - 0.22, h$lower, y + 0.22, lwd = 1.6, col = ink2)
  segments(h$upper, y - 0.22, h$upper, y + 0.22, lwd = 1.6, col = ink)
  points(h$score, y, pch = 16, cex = 1.1, col = ds_cols[h$dataset])
  cal <- (h$score - h$lower) / (h$upper - h$lower)
  mtext(sprintf("%.2f", cal), side = 4, at = y, las = 1, line = 0.3,
        cex = cs, col = ink2)
  mtext("calib.", side = 4, at = n + 0.8, las = 1, line = 0.1, cex = cs,
        col = ink2)
  title(main = "Headroom and model score", adj = 0, col.main = ink,
        cex.main = 1, font.main = 2)
  panel_label("b")
}

# (c) Linear minus training mean by reliability stratum -------------------------
panel_c <- function() {
  d <- fr_cmp[fr_cmp$method_a == "linear" &
                fr_cmp$metric %in% c("pearson_delta", "pearson_delta_top") &
                fr_cmp$condition %in% c("IFNg", "coculture"), ]
  keys <- expand.grid(stratum = c("all", "reliable", "unreliable"),
                      condition = c("IFNg", "coculture"),
                      stringsAsFactors = FALSE)
  n <- nrow(keys)
  y <- rev(seq_len(n))
  lab <- vapply(seq_len(n), function(i) {
    r <- d[d$condition == keys$condition[i] & d$stratum == keys$stratum[i], ]
    sprintf("%s\n%s (%d)", sub("Frangieh ", "", ds_labels[keys$condition[i]]),
            keys$stratum[i], if (nrow(r)) r$n_distinct[1] else 0L)
  }, character(1))
  rng <- range(c(d$ci_low, d$ci_high, 0), na.rm = TRUE)
  rng <- rng + c(-1, 1) * 0.08 * diff(rng)
  par(mar = c(3.6, 7.6, 2.6, 0.8))
  plot(NA, xlim = rng, ylim = c(0.5, n + 0.5), axes = FALSE,
       xlab = "Linear minus mean (95% CI)", ylab = "")
  xt <- pretty(rng, 5)
  grid_x(xt)
  abline(v = 0, col = ink2, lwd = 1.2)
  axis(1, at = xt, labels = tick_label(xt), lwd = 0, lwd.ticks = 1,
       col.ticks = grid_col)
  axis(2, at = y, labels = lab, las = 1, lwd = 0, cex.axis = cs, padj = 0.5)
  off <- c(pearson_delta = 0.15, pearson_delta_top = -0.15)
  for (i in seq_len(n)) {
    for (mt in names(off)) {
      r <- d[d$condition == keys$condition[i] & d$stratum == keys$stratum[i] &
               d$metric == mt, ]
      if (nrow(r) != 1 || is.na(r$diff)) next
      col <- ds_cols[keys$condition[i]]
      yy <- y[i] + off[[mt]]
      segments(r$ci_low, yy, r$ci_high, yy, lwd = 1.6, col = col)
      points(r$diff, yy, pch = if (mt == "pearson_delta") 16 else 21,
             bg = surface, col = col, cex = 1.05, lwd = 1.3)
    }
  }
  abline(h = c(3.5), col = grid_col, lwd = 1)
  legend("bottomright", legend = c("all genes", "top 20 genes"),
         pch = c(16, 21), pt.bg = surface, col = ink2, bty = "n", cex = cs,
         text.col = ink2, title = "Pearson\u0394", title.col = ink2)
  title(main = "Linear baseline vs training mean", adj = 0, col.main = ink,
        cex.main = 1, font.main = 2)
  panel_label("c")
}

# (d) PDS of the training-mean prediction vs effect size ---------------------------
pds_vs_effect <- function() {
  f <- fr_res[fr_res$method == "mean" & fr_res$metric == "pds", ]
  f <- aggregate(value ~ condition + original, data = f, FUN = mean)
  names(f) <- c("dataset", "perturbation", "pds")
  n <- nor_res[nor_res$method == "mean" & nor_res$metric == "pds" &
                 nor_res$task == "unseen_single", ]
  n <- aggregate(value ~ perturbation, data = n, FUN = mean)
  n <- data.frame(dataset = "norman", perturbation = n$perturbation,
                  pds = n$value)
  p <- rbind(f, n)
  merge(p, sig[c("dataset", "perturbation", "effect_size", "reliability")],
        by = c("dataset", "perturbation"))
}

panel_d <- function() {
  p <- pds_vs_effect()
  p <- p[is.finite(p$effect_size) & p$effect_size > 0, ]
  xr <- range(p$effect_size)
  par(mar = c(3.6, 3.8, 2.6, 0.8))
  plot(NA, xlim = xr, ylim = c(0, 1), log = "x", axes = FALSE,
       xlab = "Effect size (RMS change per gene)",
       ylab = "PDS of training-mean prediction")
  xt <- c(0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3)
  xt <- xt[xt >= xr[1] & xt <= xr[2]]
  grid_x(xt); grid_y(seq(0, 1, 0.25))
  axis(1, at = xt, labels = tick_label(xt), lwd = 0, lwd.ticks = 1,
       col.ticks = grid_col)
  axis(2, at = seq(0, 1, 0.25), labels = tick_label(seq(0, 1, 0.25)),
       las = 1, lwd = 0, lwd.ticks = 1, col.ticks = grid_col)
  abline(h = 0.5, lty = 3, col = ink2)
  text(xr[2], 0.5, "chance level", adj = c(1, 1.4), cex = cs,
       col = ink2)
  for (d in c("monoculture", "IFNg", "coculture", "norman")) {
    s <- p[p$dataset == d, ]
    points(s$effect_size, s$pds, pch = 16, cex = 0.5,
           col = alpha(ds_cols[d], 0.6))
    # smooth trend: lowess on the log effect-size scale (binned means were
    # jagged for Norman, which has only ~100 single perturbations)
    ok <- is.finite(s$effect_size) & is.finite(s$pds) & s$effect_size > 0
    if (sum(ok) >= 20) {
      fit <- stats::lowess(log(s$effect_size[ok]), s$pds[ok], f = 0.5)
      lines(exp(fit$x), pmin(pmax(fit$y, 0), 1), col = ds_cols[d], lwd = 1.8)
    }
  }
  legend("topright", legend = ds_labels, lwd = 2, col = ds_cols, bty = "n",
         cex = cs, text.col = ink2, bg = surface)
  title(main = "PDS chance level vs effect size", adj = 0,
        col.main = ink, cex.main = 1, font.main = 2)
  panel_label("d")
}

draw_fig4 <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(2, 2), oma = c(0, 0, 0, 0), bg = surface,
      col.axis = ink2, col.lab = ink2, fg = ink2, mgp = c(2.2, 0.55, 0),
      tcl = -0.3)
  par(cex = 1)   # mfrow would otherwise shrink all text to 83%
  panel_a(); panel_b(); panel_c(); panel_d()
}

# Drawn at final print size (7.5 x 7 inches, 9-point base font), so the text
# sizes are the sizes in the published figure; PNG at 400 dpi (3000 px wide).
png(file.path(out_dir, "fig4_signal_strata.png"), width = 3000, height = 2800,
    res = 400, pointsize = 9)
draw_fig4()
dev.off()
pdf_ok <- tryCatch({
  grDevices::cairo_pdf(file.path(out_dir, "fig4_signal_strata.pdf"),
                       width = 7.5, height = 7, pointsize = 9)
  draw_fig4()
  dev.off()
  TRUE
}, error = function(e) FALSE)
if (!pdf_ok) message("PDF not written (cairo_pdf unavailable); PNG is fine.")
message("Figure written to ", file.path(out_dir, "fig4_signal_strata.png"))
