#' Benchmark RTLbase training and update routines
#'
#' @description
#' Runs the baseline classifier and downstream updates on synthetic datasets of
#' varying sizes to quantify speedups from the new vectorized and parallel
#' execution paths.
#'
#' @param sizes Named integer vector of dataset sizes to benchmark.
#' @param features Integer vector of feature counts aligned to `sizes` (recycled
#'   when length-one).
#' @param repetitions Number of repeated runs per dataset size and mode.
#' @param parallel_cores Optional integer overriding detected cores for parallel
#'   runs.
#' @param margin Margin value forwarded to [alg4_BiasUpdate()].
#' @param seed Integer seed for reproducible synthetic data generation.
#' 
#' @return A `data.table` summarizing elapsed seconds per dataset and mode.
#' @export
rtl_benchmark <- function(sizes = c(small = 250, medium = 1500, large = 3200),
                          features = c(6, 24, 64),
                          repetitions = 2,
                          parallel_cores = NULL,
                          margin = 0.5,
                          seed = 99){

  seed <- normalize_seed(seed)

  required_pkgs <- c("e1071", "caret")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop(paste("Install missing packages before benchmarking:", paste(missing_pkgs, collapse = ", ")))
  }
  
  # Check optional packages that may be needed for specific algorithms
  optional_pkgs <- c("GSE", "gradDescent")
  missing_optional <- optional_pkgs[!vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_optional)) {
    warning(paste("Optional packages not available (some features may be limited):", 
                  paste(missing_optional, collapse = ", ")))
  }

  if(length(features) == 1) features <- rep(features, length(sizes))
  if(length(features) != length(sizes)) stop("Length of `features` must match `sizes`.")

  cores_used <- if (is.null(parallel_cores)) max(1L, parallel::detectCores() - 1L) else parallel_cores
  wide_threshold <- max(50, max(features))

  synth_dataset <- function(n, p){
    data.frame(matrix(rnorm(n * p), ncol = p, dimnames = list(NULL, paste0("feat", seq_len(p)))))
  }

  run_pipeline <- function(n, p, mode_label){
    TrainXls <- list(src1 = synth_dataset(n, p), src2 = synth_dataset(n, p))
    TrainYls <- lapply(TrainXls, function(df) factor(sample(c("Neg","Pos"), nrow(df), replace = TRUE), levels = c("Neg","Pos")))
    TestXls <- list(task1 = synth_dataset(n, p))
    TestYls <- list(task1 = factor(sample(c("Neg","Pos"), nrow(TestXls[[1]]), replace = TRUE), levels = c("Neg","Pos")))

    elapsed <- system.time({
      alg1 <- alg1_baselineClass(TrainXls = TrainXls,
                                 TrainYls = TrainYls,
                                 TestXls = TestXls,
                                 TestYls = TestYls,
                                 K_forCrossV = 2,
                                 svmGamma = 0.1,
                                 svmCost = 1,
                                 prnt2scr = FALSE,
                                 X_cols2Keep = seq_len(p),
                                 transX = FALSE,
                                 sampleRed = FALSE,
                                 doParalellSVM = FALSE,
                                 datatyp = "FC",
                                 use_parallel = identical(mode_label, "parallel"),
                                 parallel_cores = cores_used,
                                 wide_data_threshold = wide_threshold)

      alg2 <- alg2_rob_meanNCov(alg1$baselineSVM, print2screen = FALSE)
      alg3 <- alg3_shiftComp(source_list = TrainXls,
                             task_list = TestXls,
                             alg1_result = alg1,
                             alg2_result = alg2,
                             print2screen = FALSE,
                             save2file = FALSE,
                             maximumLag = 2,
                             ImpFeats = NA,
                             ADM = FALSE,
                             datatyp = "FC",
                             useAbsCor = TRUE,
                             medianMediansBL = FALSE,
                             CoreClassifier = "LinSVM",
                             use_parallel = identical(mode_label, "parallel"),
                             parallel_cores = cores_used,
                             wide_data_threshold = wide_threshold)

      alg4 <- alg4_BiasUpdate(task_list = TestXls,
                              alg1_result = alg1,
                              alg2_result = alg2,
                              alg3_result = alg3,
                              goodColumns = "",
                              save2file = FALSE,
                              Marg = margin,
                              alg4MinFx = "gd",
                              ADM = FALSE,
                              useMedian = TRUE,
                              ZnormMappingBL = FALSE,
                              datatyp = "FC",
                              RCSmodeBL = FALSE,
                              RCSfreqSet = c(0,0),
                              CoreClassifier = "LinSVM",
                              use_parallel = identical(mode_label, "parallel"),
                              parallel_cores = cores_used,
                              wide_data_threshold = wide_threshold)
      alg4
    })["elapsed"]

    data.table::data.table(mode = mode_label, elapsed = as.numeric(elapsed))
  }

  run_benchmark <- function(){
    results <- data.table::rbindlist(lapply(seq_along(sizes), function(idx){
    n <- sizes[[idx]]
    p <- features[[idx]]
    dataset_name <- names(sizes)[idx]

    data.table::rbindlist(lapply(seq_len(repetitions), function(rep_id){
      if (!is.null(seed)) set.seed(seed + idx + rep_id)
      seq_res <- run_pipeline(n, p, "sequential")
      seq_res[, `:=`(dataset = dataset_name, n = n, p = p, repetition = rep_id, cores = 1L)]

      par_res <- run_pipeline(n, p, "parallel")
      par_res[, `:=`(dataset = dataset_name, n = n, p = p, repetition = rep_id, cores = cores_used)]

      data.table::rbindlist(list(seq_res, par_res))
    }))
  }))
  }

  summarise_results <- function(dt) {
    dt[, .(dataset, n, p, repetition, mode, cores, elapsed)]
  }

  if (!is.null(seed)) {
    return(scoped_seed(seed, summarise_results(run_benchmark())))
  }

  summarise_results(run_benchmark())
}
