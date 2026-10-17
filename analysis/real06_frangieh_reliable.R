# Frangieh et al. 2021: evaluation restricted to reliably measured perturbations
# ---------------------------------------------------------------------------
# real05 showed that most Frangieh perturbations have split-half reliability
# near 0.1, i.e. their measured change is mostly sampling noise, and that all
# baselines then score alike. Here each condition is re-evaluated in three
# strata, defined by pe_signal() reliability BEFORE any model is fitted:
#   all, reliable (reliability >= min_rel), unreliable (< min_rel).
# Reliable perturbations are few (~2-13 per condition), so 20 random splits
# are pooled; each held-out instance is labelled "split<k>:<perturbation>".
#
# Needs the cache built by real05 (analysis/data/frangieh_cache.rds).
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real06_frangieh_reliable.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
devtools::load_all(pkg, quiet = TRUE)

data_dir <- file.path(pkg, "analysis", "data")
out_dir <- file.path(pkg, "analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(data_dir, "frangieh_cache.rds")
if (!file.exists(cache)) stop("Run real05_frangieh_baselines.R first.")

# Settings (min_cells must match real05)
min_cells <- 30
min_rel <- 0.5        # reliability threshold, fixed before evaluation
seeds <- 1:20         # random splits pooled per condition
conditions <- c("monoculture", "IFNg", "coculture")

fr <- readRDS(cache)
make_condition <- function(cn) {
  cells <- fr$condition == cn
  pert <- fr$perturbation[cells]
  n <- table(pert)
  ok <- names(n)[n >= min_cells | names(n) == "control"]
  keep <- which(cells)[pert %in% ok]
  pe_data(fr$expr[, keep], fr$perturbation[keep], control = "control")
}

# Pool per-split tables, labelling perturbations by split
pool <- function(tabs) {
  settings <- attr(tabs[[1]], "settings")
  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  attr(out, "settings") <- settings
  out
}

all_res <- list()
all_cmp <- list()
all_cal <- list()
all_per <- list()

for (cn in conditions) {
  cat("\n==", cn, "==\n")
  x <- make_condition(cn)
  sig <- pe_signal(x, seed = 1)
  per <- sig$per_perturbation
  reliable <- per$perturbation[!is.na(per$reliability) &
                                 per$reliability >= min_rel]
  cat("Reliable perturbations (reliability >=", min_rel, "):",
      length(reliable), "of", nrow(per), "\n  ",
      paste(reliable, collapse = ", "), "\n")

  res_l <- list()
  bnd_l <- list()
  for (sd in seeds) {
    s <- pe_split(x, type = "unseen_single", test_frac = 0.25, seed = sd)
    r <- pe_evaluate(list(baseline_control(x, s), baseline_mean(x, s),
                          baseline_linear(x, s)), x)
    b <- pe_bounds(x, s, seed = sd)
    r$original <- r$perturbation
    b$original <- b$perturbation
    r$perturbation <- paste0("split", sd, ":", r$perturbation)
    b$perturbation <- paste0("split", sd, ":", b$perturbation)
    res_l[[length(res_l) + 1]] <- r
    bnd_l[[length(bnd_l) + 1]] <- b
    cat(sprintf("\r  split %d/%d", sd, length(seeds)))
  }
  cat("\n")
  res <- pool(res_l)
  bnd <- pool(bnd_l)

  strata <- list(all = unique(res$perturbation),
                 reliable = unique(res$perturbation[res$original %in% reliable]),
                 unreliable = unique(res$perturbation[!res$original %in% reliable]))
  for (st in names(strata)) {
    ids <- strata[[st]]
    if (length(ids) < 2) {
      cat("  stratum", st, ": fewer than two held-out instances, skipped\n")
      next
    }
    r_st <- pe_subset(res, ids)
    b_st <- pe_subset(bnd, ids)
    n_distinct <- length(unique(r_st$original))
    cmp <- pe_compare(r_st[c("method", "perturbation", "metric", "value")],
                      reference = "mean", seed = 1)
    cal <- pe_calibrate(r_st, b_st)
    tag <- data.frame(condition = cn, stratum = st, n_distinct = n_distinct)
    all_cmp[[length(all_cmp) + 1]] <- cbind(tag, cmp)
    all_cal[[length(all_cal) + 1]] <- cbind(tag, cal)
  }

  # Per reliable perturbation: linear vs mean, averaged over the splits in
  # which it was held out
  if (length(reliable) > 0) {
    d <- res[res$original %in% reliable &
               res$metric %in% c("pearson_delta", "pearson_delta_top") &
               res$method %in% c("mean", "linear"), ]
    agg <- aggregate(value ~ original + metric + method, data = d, FUN = mean)
    wide <- reshape(agg, idvar = c("original", "metric"), timevar = "method",
                    direction = "wide")
    names(wide) <- sub("^value\\.", "", names(wide))
    held <- table(res$original[res$method == "mean" &
                                 res$metric == "pearson_delta" &
                                 res$original %in% reliable])
    wide$times_held_out <- as.integer(held[wide$original])
    wide$reliability <- per$reliability[match(wide$original, per$perturbation)]
    all_per[[length(all_per) + 1]] <- cbind(condition = cn, wide)
  }
  res$condition <- cn
  all_res[[length(all_res) + 1]] <- res
}

cmp <- do.call(rbind, all_cmp)
cal <- do.call(rbind, all_cal)
per_rel <- do.call(rbind, all_per)
write.csv(do.call(rbind, all_res),
          file.path(out_dir, "frangieh_reliable_per_perturbation.csv"),
          row.names = FALSE)
write.csv(cmp, file.path(out_dir, "frangieh_reliable_compare_vs_mean.csv"),
          row.names = FALSE)
write.csv(cal, file.path(out_dir, "frangieh_reliable_calibrated.csv"),
          row.names = FALSE)
write.csv(per_rel, file.path(out_dir, "frangieh_reliable_by_perturbation.csv"),
          row.names = FALSE)

# Printouts -------------------------------------------------------------------------
focus <- c("pearson_delta", "pearson_delta_top", "mse_top", "pds")

cat("\n== Linear vs training mean, by stratum (diff = linear - mean; 20 splits pooled) ==\n")
sc <- cmp[cmp$method_a == "linear" & cmp$metric %in% focus,
          c("condition", "stratum", "metric", "n", "n_distinct", "mean_b",
            "mean_a", "diff", "ci_low", "ci_high", "p_adj", "favours")]
names(sc)[names(sc) == "mean_b"] <- "mean_score"
names(sc)[names(sc) == "mean_a"] <- "linear_score"
num <- c("mean_score", "linear_score", "diff", "ci_low", "ci_high")
sc[num] <- round(sc[num], 3)
sc$p_adj <- signif(sc$p_adj, 2)
print(sc[order(sc$condition, sc$metric, sc$stratum), ], row.names = FALSE)

cat("\n== Calibrated score of the linear model, by stratum ==\n")
cl <- cal[cal$method == "linear" & cal$metric %in% focus,
          c("condition", "stratum", "metric", "n", "score", "lower", "upper",
            "calibrated")]
cl[c("score", "lower", "upper", "calibrated")] <-
  round(cl[c("score", "lower", "upper", "calibrated")], 3)
print(cl[order(cl$condition, cl$metric, cl$stratum), ], row.names = FALSE)

cat("\n== Reliable perturbations: mean vs linear (averaged over splits where held out) ==\n")
pr <- per_rel
pr[c("mean", "linear", "reliability")] <-
  round(pr[c("mean", "linear", "reliability")], 3)
print(pr[order(pr$condition, pr$metric, -pr$reliability),
         c("condition", "metric", "original", "reliability", "times_held_out",
           "mean", "linear")], row.names = FALSE)

cat("\nResults written to", out_dir, "\n")
