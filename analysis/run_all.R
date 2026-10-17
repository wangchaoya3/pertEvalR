# Re-run every analysis of the manuscript with the installed package version
# ---------------------------------------------------------------------------
# Order: simulation (Figure 1), Norman (Figure 3 data), Figure 3, Frangieh,
# reliability strata, Figure 4, meta-metrics. Needs the data caches in
# analysis/data/ (built on first run by real02 and real05 from the scPerturb
# files; see real01 and real04 for the downloads).
# Takes roughly 30-60 minutes.
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "run_all.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
scripts <- c("fig1_systematic_variation.R", "real02_norman_baselines.R",
             "real03_calibration_figure.R", "real05_frangieh_baselines.R",
             "real06_frangieh_reliable.R", "real07_fig4_signal.R",
             "real08_meta_metrics.R")
log_file <- file.path(pkg, "analysis", "output", "run_all_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
cat("pertEvalR", unname(read.dcf(file.path(pkg, "DESCRIPTION"))[, "Version"]),
    "run started", format(Sys.time()), "\n", file = log_file)
for (f in scripts) {
  message("\n==== ", f, " ====")
  t0 <- Sys.time()
  source(file.path(pkg, "analysis", f), echo = FALSE)
  cat(f, "finished in", round(difftime(Sys.time(), t0, units = "mins"), 1),
      "min\n", file = log_file, append = TRUE)
}
cat("All analyses finished", format(Sys.time()), "\n", file = log_file,
    append = TRUE)
message("Done. Log: ", log_file)
