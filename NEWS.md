# pertEvalR 0.0.0.9005

* New `pe_signal()`: dataset diagnostics before any model is evaluated.
  Per perturbation: cosine similarity of its change to the average change
  (systematic variation), effect size, effective number of changed genes
  (participation ratio) and split-half reliability (Spearman-Brown). Dataset
  summary: mean / median cosine and the systematic fraction of the total
  squared change.
* Internal: `pe_bounds()` and `pe_signal()` share one cell-halving helper.

# pertEvalR 0.0.0.9004

* `pe_bounds()` lower bound: the better uninformative baseline is now chosen
  once per metric (by its mean over perturbations) instead of separately for
  each perturbation, which biased the bound upwards.
* `pe_bounds()` upper bound: the split-half replicate score is corrected to
  full sample size (MSE divided by 4; Spearman-Brown ceiling
  `sqrt(2r / (1 + r))` for correlation metrics). The uncorrected score is kept
  in the new column `upper_raw`.

# pertEvalR 0.0.0.9003

* New `pe_bounds()`: per-perturbation lower bound (best of the no-change and
  training-mean baselines) and upper bound (split-half technical replicate)
  for every metric.
* New `pe_calibrate()`: places method scores between those bounds
  (0 = uninformative, 1 = technical replicate).
* `pe_evaluate()` results now record their `top_n` and `distance` settings in
  a `"settings"` attribute, so mismatched calibrations are caught.

# pertEvalR 0.0.0.9002

* New `pe_compare()`: paired comparison of methods per metric, with the mean
  per-perturbation difference, a percentile bootstrap confidence interval,
  a paired Wilcoxon signed-rank test with multiple-testing correction, and a
  `favours` verdict that respects each metric's direction.

# pertEvalR 0.0.0.9001

* New metrics via `pe_evaluate()`: `mse`, `mse_top`, `pearson_delta`,
  `pearson_delta_top`, `pearson_delta_pert` (centroid reference) and `pds`
  (perturbation discrimination score); `pe_summary()` averages them.
* New baselines `baseline_additive()` (combinations) and `baseline_linear()`
  (Ahlmann-Eltze et al. 2025 linear model).
* New `pe_prediction()` wraps predictions from any model for evaluation.
* `pe_simulate()` now gives each single perturbation a target gene
  (`x$sim$targets`).
* Tests for the SingleCellExperiment and Seurat converters (skipped when
  those packages are not installed).

# pertEvalR 0.0.0.9000

* Initial skeleton: data layer (`pe_data()`, `as_pe_data()`, `pe_pseudobulk()`,
  `pe_split()`, `pe_simulate()`) and two reference baselines
  (`baseline_control()`, `baseline_mean()`).
