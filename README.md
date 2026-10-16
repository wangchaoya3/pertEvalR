# pertEvalR

> Calibrated evaluation of single-cell perturbation response predictions in R.
> **Status: early skeleton (v0.0.0.9000).** The name is a placeholder; check
> CRAN, Bioconductor and GitHub for clashes before the first release.

Models that predict how cells respond to perturbations ("virtual cell" and
in-silico perturbation models) are increasingly compared using metrics that
uninformative baselines can score well on. `pertEvalR` aims to give R and
Bioconductor users one framework for:

- standard, leakage-free **splits** of perturbations;
- **reference baselines** that any useful model must beat;
- **calibrated metrics** (scores placed between an uninformative lower bound
  and a technical-replicate upper bound);
- **signal diagnostics** (how much systematic variation a dataset carries, and
  whether a metric can tell models apart on it);
- **statistical comparison** of models (confidence intervals and paired tests
  across perturbations), which existing evaluation tools largely lack.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YOUR-GITHUB/pertEvalR")
```

## Quick start

```r
library(pertEvalR)

x <- pe_simulate(seed = 1)          # or pe_data(expr, labels) / as_pe_data(sce, ...)
s <- pe_split(x, type = "unseen_single", seed = 42)

truth  <- pe_truth(x, s)            # observed pseudobulk of held-out perturbations
p_ctrl <- baseline_control(x, s)    # "no change"
p_mean <- baseline_mean(x, s)       # average training response

mean((p_ctrl$pred - truth)^2)
mean((p_mean$pred - truth)^2)
```

From a `SingleCellExperiment` or a Seurat v5 object:

```r
x <- as_pe_data(sce, perturbation_col = "target_gene", control = "non-targeting")
x <- as_pe_data(seu, perturbation_col = "target_gene", control = "non-targeting")
```

Predictions from your own model only need to be a genes x held-out
perturbations matrix laid out like `pe_truth(x, s)`.

## What is implemented

| Layer | Function | Status |
|---|---|---|
| Data | `pe_data()`, `as_pe_data()` (matrix, SCE, Seurat v5) | done |
| Data | `pe_pseudobulk()`, `pe_truth()`, `pe_perturbations()` | done |
| Data | `pe_split()` — `unseen_single`, `unseen_combo` | done |
| Data | `pe_simulate()` — toy data with known effects | done |
| Baselines | `baseline_control()`, `baseline_mean()` | done |

## Roadmap

| Layer | Planned | Notes |
|---|---|---|
| Data | `unseen_cell_line` split; h5ad import via `zellkonverter` | |
| Baselines | `baseline_additive()`, `baseline_linear()` | linear model after Ahlmann-Eltze et al. (Nat Methods 2025) |
| Baselines | `baseline_linear_prior()` | external gene embeddings, e.g. STRING |
| Metrics | MSE/MAE (+ DEG-weighted), PearsonΔ vs control and vs perturbed centroid, PDS, DE overlap/precision/direction, energy distance, MMD | centroid reference follows Systema (Nat Biotech 2025) |
| Calibration | `pe_bounds()` — split-half upper bound, baseline lower bound | after the SBB principles (bioRxiv 2026) |
| Diagnostics | `pe_signal()` — systematic variation, effective gene number | |
| Statistics | `pe_compare()` — paired bootstrap / Wilcoxon across perturbations | main point of difference from existing tools |
| Reporting | `pe_report()` HTML; `pe_parity_celleval()` | parity with Arc's cell-eval |
| Docs | vignette on Norman 2019 and Replogle 2022 | |

## Development

```r
devtools::document()   # generates man/ and refreshes NAMESPACE
devtools::test()
devtools::check()
```

Commit the generated `man/` folder; the GitHub Actions workflow in
`.github/workflows/` runs `R CMD check` on Linux, macOS and Windows.

## License

MIT
