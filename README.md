# RTLbase
## The R Transfer Learning (RTL) Classification Framework - Baseline Implementation v1.0.1

<a href="https://zenodo.org/badge/latestdoi/160988086"><img src="https://zenodo.org/badge/160988086.svg" alt="DOI"></a>


-----------

Author/Maintainer: Eisa Mahyari

Twitter: @eisamahyari

URL: https://github.com/eisascience/RTLbase

-----------
Current Development: RTLbase v1.0.1, a sub-package of RTL v1.0.0

A transductive transfer learning (TL) classification framework, that adapts generalized linear classification hyperplane obtained on training data to the previously unseen testing data to identify distinct classes; specifically, subsets of cells defined by a phenotype or a signature from flow-based and transcriptomic single-cell assays.

This framework was developed as part of a PhD dissertation research, focused on developing robust and reproducible methods to identify rare cellular subsets from high-dimensional single-cell data (e.g., CyTOF or scRNASeq). 

-------------

Based on a workshop manuscript proposed by Lee et al. (2011) for automated classification of flowcytometry data

-------------

## Quickstart

1. Install the development dependencies:

```r
install.packages(c("devtools", "roxygen2", "knitr", "rmarkdown"))
```

2. Install RTLbase from a local checkout or directly from GitHub:

```r
# from a local clone
devtools::install_local(".")

# or install the latest main branch
remotes::install_github("eisascience/RTLbase")
```

3. Load the package and verify the main functions are available:

```r
library(RTLbase)
ls("package:RTLbase")
```

## Minimal example

The snippet below generates a small synthetic dataset, fits the baseline RTL
models, and produces summary statistics for the adapted hyperplanes. The code is
kept lightweight so you can paste it into an R console and see the workflow.

```r
library(RTLbase)

# simulate three train/test splits with two features each
source("R/simData.R")
demo <- DemoTestTrainList2D(N = 3)

train_ls <- lapply(demo, function(x) list(X.train = x$TrainX, Y.train = x$TrainY))
test_ls  <- lapply(demo, function(x) list(X.test = x$TestX,  Y.test = x$TestY))

alg1 <- alg1_baselineClass(
  TrainXls = lapply(train_ls, `[[`, "X.train"),
  TrainYls = lapply(train_ls, `[[`, "Y.train"),
  TestXls  = lapply(test_ls, `[[`, "X.test"),
  TestYls  = lapply(test_ls, `[[`, "Y.test"),
  K_forCrossV = 3,
  svmGamma = 0.02,
  svmCost = 0.5,
  verbose = TRUE,
  X_cols2Keep = 1:2,
  transX = FALSE,
  sampleRed = FALSE,
  doParalellSVM = FALSE
)

alg2 <- alg2_rob_meanNCov(alg1$baselineSVM)
alg3 <- alg3_shiftComp(
  source_list = lapply(train_ls, `[[`, "X.train"),
  task_list = lapply(test_ls, `[[`, "X.test"),
  alg1_result = alg1,
  alg2_result = alg2,
  verbose = FALSE,
  save2file = FALSE,
  maximumLag = 0,
  ImpFeats = NA,
  CoreClassifier = "LinSVM"
)

alg4 <- alg4_BiasUpdate(
  task_list = lapply(test_ls, `[[`, "X.test"),
  alg1_result = alg1,
  alg2_result = alg2,
  alg3_result = alg3,
  goodColumns = NA,
  save2file = FALSE,
  Marg = 1,
  alg4MinFx = "gd",
  CoreClassifier = "LinSVM"
)

alg6 <- alg6_NormalVectorUpdate(
  task_list = lapply(test_ls, `[[`, "X.test"),
  alg1_result = alg1,
  alg2_result = alg2,
  alg3_result = alg3,
  alg4_result = alg4,
  X_feat_cols = NA,
  Marg = 1,
  save2file = FALSE,
  ADM = FALSE,
  datatyp = "FC",
  CoreClassifier = "LinSVM"
)

viz <- FinalViz(
  TrainTestSet.ls = test_ls,
  alg1_result = alg1,
  alg2_result = alg2,
  alg3_result = alg3,
  alg4_result = alg4,
  alg6_result = alg6,
  datatyp = "FC",
  ADM = FALSE
)

viz$comStats.sub
```

## Performance notes

- Start with the smallest number of simulated splits or down-sample (`sampleRed`
  in `alg1_baselineClass`) when prototyping new pipelines.
- Use `doParalellSVM = TRUE` only after verifying memory availability; large
  flow cytometry datasets can easily exhaust RAM.
- The `maximumLag` parameter in `alg3_shiftComp` controls the lag search window;
  keeping it near `0` speeds up experimentation while larger windows improve
  alignment on challenging datasets.
- Set `save2file = TRUE` when running long experiments so that intermediate
  plots and cross-validation summaries are persisted for later review.
- New `use_parallel` and `parallel_cores` arguments on the classifier and
  shift-compensation routines batch work across datasets using `parallel::mclapply`;
  set `wide_data_threshold` to coerce wide feature matrices through
  `data.table` to reduce copy overhead.
- Run `rtl_benchmark()` to compare sequential versus parallel execution on
  synthetic small/medium/large workloads:

```r
benchmark_results <- rtl_benchmark()
print(benchmark_results)
```
