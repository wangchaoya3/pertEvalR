# Assemble the Zenodo / PeerJ release bundle for the current version
# ---------------------------------------------------------------------------
# Run AFTER: (1) run_all.R has regenerated all outputs with this version,
# (2) all changes are committed, (3) the version tag (e.g. v0.1.0) exists.
# Writes to Desktop/pertEvalR_zenodo_<version>/:
#   pertEvalR_<version>.tar.gz         R source package (R CMD build)
#   pertEvalR-<version>-source.zip     every file tracked by git at this commit
#                                      (package + analysis scripts)
#   pertEvalR-<version>-results.zip    all result tables (CSV) behind the figures
#   figures/Figure 1-4 .png/.pdf       manuscript figures, PeerJ naming
#   sessionInfo.txt                    R and package versions
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "make_release.R"))

pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
ver <- unname(read.dcf(file.path(pkg, "DESCRIPTION"))[, "Version"])
tag <- paste0("v", ver)
if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")

# 0. Checks ------------------------------------------------------------------
st <- gert::git_status(repo = pkg)
if (nrow(st) > 0) {
  print(st)
  stop("Uncommitted changes: commit them first so the archive equals the tag.")
}
tags <- gert::git_tag_list(repo = pkg)
if (!tag %in% tags$name) stop("Tag ", tag, " not found; create it first.")
head_commit <- gert::git_commit_id("HEAD", repo = pkg)
tag_commit <- gert::git_commit_id(tag, repo = pkg)
if (!identical(head_commit, tag_commit)) {
  stop("HEAD is not at ", tag, "; check out the tag or move it.")
}

out <- file.path(dirname(pkg), paste0("pertEvalR_zenodo_", ver))
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
res_dir <- file.path(pkg, "analysis", "output")

# 1. R source package ----------------------------------------------------------
R <- file.path(R.home("bin"), "R")
old <- setwd(dirname(pkg))
system2(R, c("CMD", "build", "pertEvalR"))
tarball <- paste0("pertEvalR_", ver, ".tar.gz")
file.copy(tarball, file.path(out, tarball), overwrite = TRUE)
setwd(old)

# 2. Full source snapshot (tracked files only) ---------------------------------
tracked <- gert::git_ls(repo = pkg)$path
zip::zip(file.path(out, paste0("pertEvalR-", ver, "-source.zip")),
         files = tracked, root = pkg)

# 3. Figures with PeerJ names -----------------------------------------------------
figs <- list(
  "Figure 1" = file.path(res_dir, "fig1_systematic_variation"),
  "Figure 2" = file.path(pkg, "analysis", "figure2", "fig2_workflow"),
  "Figure 3" = file.path(res_dir, "fig3_calibrated_scores"),
  "Figure 4" = file.path(res_dir, "fig4_signal_strata")
)
png_size <- function(f) {
  b <- readBin(f, "raw", 24)
  w <- sum(as.integer(b[17:20]) * 256^(3:0))
  h <- sum(as.integer(b[21:24]) * 256^(3:0))
  c(width = w, height = h)
}
for (nm in names(figs)) {
  for (ext in c("png", "pdf")) {
    src <- paste0(figs[[nm]], ".", ext)
    if (file.exists(src)) {
      file.copy(src, file.path(out, "figures", paste0(nm, ".", ext)),
                overwrite = TRUE)
    } else {
      warning("Missing ", src)
    }
  }
  p <- paste0(figs[[nm]], ".png")
  if (file.exists(p)) {
    sz <- png_size(p)
    cat(sprintf("%s: %d x %d px%s\n", nm, sz[1], sz[2],
                if (max(sz) < 900) "  <-- below PeerJ minimum of 900 px" else ""))
  }
}

# 4. Result tables ------------------------------------------------------------------
csvs <- list.files(res_dir, pattern = "\\.csv$")
zip::zip(file.path(out, paste0("pertEvalR-", ver, "-results.zip")),
         files = c(csvs, intersect("run_all_log.txt", list.files(res_dir))),
         root = res_dir)

# 5. Session information -------------------------------------------------------------
writeLines(c(paste("pertEvalR", ver, "commit", head_commit),
             capture.output(sessionInfo())),
           file.path(out, "sessionInfo.txt"))

cat("\nRelease bundle written to", out, "\n")
print(data.frame(file = list.files(out, recursive = TRUE),
                 MB = round(file.size(list.files(out, recursive = TRUE,
                                                 full.names = TRUE)) / 1e6, 2)))
