# Meta-metrics (dynamic range fraction, bound discrimination score) per dataset
# ---------------------------------------------------------------------------
# For every dataset/task/condition: pe_bounds() on 5 random splits (seeds 1-5,
# 25% held out), pooled, then pe_meta_metrics(). For Frangieh also restricted
# to reliable perturbations (pe_signal reliability >= 0.5).
# Needs the caches built by real02 and real05.
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real08_meta_metrics.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
devtools::load_all(pkg, quiet = TRUE)
data_dir <- file.path(pkg, "analysis", "data")
out_dir <- file.path(pkg, "analysis", "output")
seeds <- 1:5
min_cells <- 30
min_rel <- 0.5

pooled_bounds <- function(x, type) {
  tabs <- lapply(seeds, function(sd) {
    s <- pe_split(x, type = type, test_frac = 0.25, seed = sd)
    b <- pe_bounds(x, s, seed = sd)
    b$original <- b$perturbation
    b$perturbation <- paste0("split", sd, ":", b$perturbation)
    b
  })
  settings <- attr(tabs[[1]], "settings")
  out <- do.call(rbind, tabs)
  attr(out, "settings") <- settings
  out
}

rows <- list()
add <- function(dataset, subset, mm) {
  rows[[length(rows) + 1]] <<- cbind(dataset = dataset, subset = subset, mm)
}

# Norman
nor <- readRDS(file.path(data_dir, "norman_pe_data.rds"))
for (tk in c("unseen_single", "unseen_combo")) {
  b <- pooled_bounds(nor, tk)
  add(paste("Norman", tk), "all", pe_meta_metrics(b))
  cat("Norman", tk, "done\n")
}

# Frangieh, per condition
fr <- readRDS(file.path(data_dir, "frangieh_cache.rds"))
for (cn in c("monoculture", "IFNg", "coculture")) {
  cells <- fr$condition == cn
  pert <- fr$perturbation[cells]
  n <- table(pert)
  ok <- names(n)[n >= min_cells | names(n) == "control"]
  keep <- which(cells)[pert %in% ok]
  x <- pe_data(fr$expr[, keep], fr$perturbation[keep], control = "control")
  sig <- pe_signal(x, seed = 1)
  per <- sig$per_perturbation
  reliable <- per$perturbation[!is.na(per$reliability) &
                                 per$reliability >= min_rel]
  b <- pooled_bounds(x, "unseen_single")
  add(paste("Frangieh", cn), "all", pe_meta_metrics(b))
  ids <- unique(b$perturbation[b$original %in% reliable])
  if (length(ids) > 0) {
    add(paste("Frangieh", cn), "reliable", pe_meta_metrics(pe_subset(b, ids)))
  }
  cat("Frangieh", cn, "done\n")
}

mm <- do.call(rbind, rows)
rownames(mm) <- NULL
write.csv(mm, file.path(out_dir, "meta_metrics.csv"), row.names = FALSE)

cat("\n== Dynamic range fraction (DRF) and bound discrimination score (BDS) ==\n")
show <- mm
show[c("lower", "upper", "drf", "bds")] <- round(show[c("lower", "upper", "drf", "bds")], 3)
print(show[order(show$metric, show$dataset, show$subset),
           c("metric", "dataset", "subset", "n", "lower", "upper", "drf", "bds")],
      row.names = FALSE)
