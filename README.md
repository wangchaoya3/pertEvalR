# pertEvalR

> Calibrated evaluation of single-cell perturbation response predictions in R.
> **Version 0.1.0.**

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
remotes::install_github("wangchaoya3/pertEvalR")
```

## Quick start

```r
library(pertEvalR)

x <- pe_simulate(seed = 1)          # or pe_data(expr, labels) / as_pe_data(sce, ...)
pe_signal(x, seed = 1)              # how much of the change is one shared shift?
s <- pe_split(x, type = "unseen_single", seed = 42)

preds <- list(
  baseline_control(x, s),                          # "no change"
  baseline_mean(x, s),                             # average training response
  baseline_linear(x, s, targets = x$sim$targets)   # Ahlmann-Eltze linear model
)
res <- pe_evaluate(preds, x)        # one row per method x perturbation x metric
pe_summary(res)

# Is the linear model really better than the training mean?
pe_compare(res, reference = "mean", seed = 1)

# How far between "uninformative" (0) and "technical replicate" (1)?
pe_calibrate(res, pe_bounds(x, s, seed = 1))
```

The vignette walks through the whole workflow, including signal
diagnostics and reliability strata:
`remotes::install_github("wangchaoya3/pertEvalR", build_vignettes = TRUE)`,
then `vignette("pertEvalR")`.

Combinations use their own split and the additive baseline:

```r
s2 <- pe_split(x, type = "unseen_combo", seed = 42)
pe_summary(pe_evaluate(baseline_additive(x, s2), x))
```

From a `SingleCellExperiment` or a Seurat v5 object:

```r
x <- as_pe_data(sce, perturbation_col = "target_gene", control = "non-targeting")
x <- as_pe_data(seu, perturbation_col = "target_gene", control = "non-targeting")
```

Predictions from your own model: wrap a genes x held-out perturbations
matrix with `pe_prediction(pred, s, method = "my_model")` and pass it to
`pe_evaluate()`.

## What is implemented

| Layer | Function | Status |
|---|---|---|
| Data | `pe_data()`, `as_pe_data()` (matrix, SCE, Seurat v5) | done |
| Data | `pe_pseudobulk()`, `pe_truth()`, `pe_perturbations()` | done |
| Data | `pe_split()` — `unseen_single`, `unseen_combo` | done |
| Data | `pe_simulate()` — toy data with known effects and targets | done |
| Baselines | `baseline_control()`, `baseline_mean()`, `baseline_additive()`, `baseline_linear()` | done |
| Predictions | `pe_prediction()` — wrap any model's output | done |
| Metrics | `pe_evaluate()`: MSE (all / top-n genes), PearsonΔ vs control (all / top-n), PearsonΔ vs perturbed centroid, PDS (L1 / L2 / cosine); `pe_summary()` | done |
| Statistics | `pe_compare()`: paired bootstrap CI, Wilcoxon signed-rank, multiple-testing correction, direction-aware verdict | done |
| Calibration | `pe_bounds()` (best uninformative baseline / split-half replicate), `pe_calibrate()`; `pe_meta_metrics()` (dynamic range fraction, bound discrimination score) | done |
| Diagnostics | `pe_signal()`: systematic variation (cosine to average shift, systematic fraction), effect size, effective gene number, split-half reliability; `pe_subset()` to evaluate only reliable perturbations | done |

## Roadmap

| Layer | Planned | Notes |
|---|---|---|
| Data | `unseen_cell_line` split; h5ad import via `zellkonverter` | |
| Baselines | `baseline_linear_prior()` | external gene embeddings, e.g. STRING |
| Metrics | MAE, DEG-weighted variants, DE overlap/precision/direction, energy distance, MMD | |
| Diagnostics | metric discriminability (can a metric separate known-good from uninformative predictions?) | |
| Statistics | across-split variation; comparisons pooled over several splits | |
| Reporting | `pe_report()` HTML; `pe_parity_celleval()` | parity with Arc's cell-eval |
| Docs | real-data article (Norman 2019, Frangieh 2021) | |

## Development

```r
devtools::document()   # generates man/ and refreshes NAMESPACE
devtools::test()
devtools::check()
```

Commit the generated `man/` folder; the GitHub Actions workflow in
`.github/workflows/` runs `R CMD check` on Linux.

## License

MIT
