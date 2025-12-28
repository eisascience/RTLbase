# RTLbase 1.0.0

## New features

* Added comprehensive package metadata, licensing, and testthat configuration.
* Introduced snapshot-based tests for deterministic helper utilities and the baseline SVM classifier structure.
* Added continuous integration via R CMD check workflow and linting/styling configuration to support consistent development.

## Performance improvements

* Added optional parallel backends and wide-data coercion across the baseline classifier and shift compensation routines to speed up batch processing.
* Vectorized bias-count aggregation to reduce per-task runtime in bias updates.
* Introduced `rtl_benchmark()` to profile small, medium, and large synthetic datasets with sequential versus parallel execution, along with a repeatable configuration for measuring gains.
