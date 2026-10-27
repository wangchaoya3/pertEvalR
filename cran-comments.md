## Submission

This is the first CRAN submission of pertEvalR.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission.
  Possibly misspelled words in DESCRIPTION: words such as "Perturb",
  "pseudobulk" and author names are correct technical terms and names.

## Test environments

* local: Windows 11, R 4.6.1
* win-builder: R-devel
* GitHub Actions: ubuntu-latest, R release

## Notes for the reviewer

* Packages in Suggests from Bioconductor (SingleCellExperiment,
  SummarizedExperiment) are used only by the optional `as_pe_data()` methods
  for those classes; tests that need them are skipped when they are not
  installed, and the vignette does not evaluate that code.
* Functions that accept a `seed` argument set the random number generator
  temporarily and restore the user's `.Random.seed` afterwards.
