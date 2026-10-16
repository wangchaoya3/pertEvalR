# Real data, step 1: download and inspect Norman et al. 2019 (scPerturb)
# ---------------------------------------------------------------------------
# Purpose: find out how the file is laid out (which obs column holds the
# perturbation label, how controls and combinations are written, matrix
# size) BEFORE writing the analysis, and check available memory.
# Nothing here loads the full expression matrix.
#
# Data: scPerturb harmonised Perturb-seq collection, Zenodo record 10044268,
# CC BY 4.0. File: NormanWeissman2019_filtered.h5ad (about 700 MB).
#
# One-time setup (Bioconductor; no Python needed):
#   install.packages("BiocManager")
#   options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#   BiocManager::install(c("anndataR", "rhdf5"))
#
# Run:
#   source(file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR",
#                    "analysis", "real01_inspect_norman.R"))
# Then paste everything printed in the console back into the chat.

# Package folder on this computer's Desktop (works for any Windows user name).
# If your Desktop is elsewhere (e.g. OneDrive), set pkg to the full path.
pkg <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "pertEvalR")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("pertEvalR not found at ", pkg, ". Set `pkg` to the package folder.")
}
data_dir <- file.path(pkg, "analysis", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

file_name <- "NormanWeissman2019_filtered.h5ad"
h5ad <- file.path(data_dir, file_name)
url <- paste0("https://zenodo.org/records/10044268/files/", file_name,
              "?download=1")

# 1. Download (skipped if the file is already there) ---------------------------
# If this is slow or fails, download the URL above in a browser and put the
# file in analysis/data/ with exactly this name, then re-run the script.
if (!file.exists(h5ad)) {
  options(timeout = 7200)
  message("Downloading ", file_name, " (~700 MB) ...")
  download.file(url, h5ad, mode = "wb")
}
cat("File size:", round(file.size(h5ad) / 1024^2), "MB\n")

# 2. Memory on this computer ----------------------------------------------------
mem <- tryCatch(ps::ps_system_memory(), error = function(e) NULL)
if (!is.null(mem)) {
  cat("RAM total:", round(mem$total / 1024^3, 1), "GB; available:",
      round(mem$avail / 1024^3, 1), "GB\n")
}

# 3. HDF5 layout (top two levels only) -------------------------------------------
cat("\n== HDF5 layout ==\n")
layout <- rhdf5::h5ls(h5ad, recursive = 2)
print(layout[, c("group", "name", "otype", "dclass", "dim")], row.names = FALSE)

# Number of stored non-zero values in X, if X is sparse (memory estimate)
nnz <- tryCatch({
  d <- layout$dim[layout$group == "/X" & layout$name == "data"]
  as.numeric(d)
}, error = function(e) NA)
if (length(nnz) == 1 && !is.na(nnz)) {
  cat("\nNon-zero values in X:", format(nnz, big.mark = ","),
      "-> about", round(nnz * 12 / 1024^3, 1),
      "GB as an R sparse matrix\n")
}

# 4. Cell and gene metadata (lazy read: X is not loaded) ------------------------
cat("\n== Metadata ==\n")
ad <- anndataR::read_h5ad(h5ad, as = "HDF5AnnData")
obs <- ad$obs
var <- ad$var
cat("Cells:", nrow(obs), " Genes:", nrow(var), "\n\n")

cat("obs columns:\n")
print(vapply(obs, function(col) class(col)[1], character(1)))
cat("\nFirst rows of obs:\n")
print(utils::head(obs, 3))

cat("\nvar columns:", paste(names(var), collapse = ", "), "\n")
cat("First gene names:", paste(utils::head(rownames(var), 5),
                               collapse = ", "), "\n")

# 5. Candidate perturbation columns ---------------------------------------------
cat("\n== Candidate label columns ==\n")
for (col in names(obs)) {
  v <- obs[[col]]
  if (!(is.character(v) || is.factor(v))) next
  n_lev <- length(unique(v))
  if (n_lev < 2 || n_lev > 2000) next
  tab <- sort(table(v), decreasing = TRUE)
  cat("\n--", col, ":", n_lev, "distinct values. Most common:\n")
  print(utils::head(tab, 5))
  ex <- utils::head(setdiff(names(tab), names(utils::head(tab, 5))), 8)
  cat("   other examples:", paste(ex, collapse = " | "), "\n")
  for (sep in c("_", "+", ",", ";", "|")) {
    k <- sum(grepl(sep, names(tab), fixed = TRUE))
    if (k > 0) cat("   labels containing '", sep, "': ", k, "\n", sep = "")
  }
}

for (col in intersect(c("nperts", "perturbation_type", "celltype",
                        "cell_line"), names(obs))) {
  cat("\n--", col, ":\n")
  print(utils::head(sort(table(obs[[col]]), decreasing = TRUE), 6))
}

cat("\nDone. Please paste all of the output above into the chat.\n")
