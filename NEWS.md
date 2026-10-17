# pertEvalR 0.1.0

First release, accompanying the manuscript submitted to PeerJ.

* Data: `pe_data()`, `as_pe_data()` (matrix, SingleCellExperiment, Seurat v5),
  `pe_pseudobulk()`, `pe_truth()`, `pe_simulate()`.
* Diagnostics: `pe_signal()` (systematic variation, effect size, effective
  genes, split-half reliability).
* Splits: `pe_split()` (unseen single perturbations, unseen combinations).
* Baselines: `baseline_control()`, `baseline_mean()`, `baseline_additive()`,
  `baseline_linear()`; `pe_prediction()` for external models.
* Metrics: `pe_evaluate()`, `pe_summary()`.
* Statistics: `pe_compare()`.
* Calibration: `pe_bounds()`, `pe_calibrate()`, `pe_meta_metrics()`.
* Stratification: `pe_subset()`.
* Vignette covering the full workflow; 288 unit tests.

# pertEvalR 0.0.0.9008

* New `pe_meta_metrics()`: dynamic range fraction (Miller et al. 2026) and
  bound discrimination score (Vollenweider and Buehlmann 2026) for every
  metric, computed from `pe_bounds()`, to judge whether a metric can separate
  informative from uninformative predictions on a dataset.

# pertEvalR 0.0.0.9007

* New vignette "Evaluating perturbation response predictions with pertEvalR"
  covering the full workflow on simulated data
  (`vignette("pertEvalR")`).
* `pe_simulate()` gains `effect_spread`: log-normal per-perturbation
  strength, so simulated screens mix strong and barely measurable
  perturbations. The default (0) reproduces earlier simulations exactly; the
  strengths are returned in `sim$strength`.

# pertEvalR 0.0.0.9006

* New `pe_subset()`: restricts `pe_evaluate()` results or `pe_bounds()` tables
  to chosen perturbations (or to those a `pe_signal()` object marks as
  reliable) while keeping the settings that `pe_calibrate()` checks.
* `?pe_signal` explains how noise affects `effective_genes` and how cell
  numbers affect `reliability`.

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
