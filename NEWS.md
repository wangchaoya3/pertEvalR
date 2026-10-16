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
