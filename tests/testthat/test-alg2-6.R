skip_if_missing_alg_deps <- function() {
  skip_if_not_installed("GSE")
  skip_if_not_installed("gradDescent")
}

build_baseline_matrix <- function() {
  mat <- matrix(c(
    1.0, 0.2,
    1.2, 0.1,
    0.9, 0.3
  ), byrow = TRUE, ncol = 2)

  colnames(mat) <- c("w1", "b.int")
  mat
}

build_small_tasks <- function() {
  list(
    src1 = data.frame(x = seq(-1, 1, length.out = 12)),
    src2 = data.frame(x = seq(0, 2, length.out = 12))
  )
}

build_target_tasks <- function() {
  list(task1 = data.frame(x = seq(-0.5, 1.5, length.out = 12)))
}

test_that("alg2 robust mean/cov returns normalized statistics", {
  skip_if_missing_alg_deps()
  alg2_res <- alg2_rob_meanNCov(build_baseline_matrix(), print2screen = FALSE)

  expect_named(alg2_res$U_robust_norm, c("w1", "b.int"))
  expect_length(alg2_res$v_0, 1)
  expect_gt(alg2_res$w_euc_mag, 0)
  expect_false(any(is.na(unlist(alg2_res[c("U_robust", "C_robust", "U_robust_norm")]))))
})

test_that("bandwidth estimator surfaces recycling and missing data", {
  skip_if_missing_alg_deps()
  expect_warning(KBand_fx(1:3, c(1, 2)), "length")
  expect_true(is.na(KBand_fx(c(1, NA, 3), c(1, 1, 1))))
})

test_that("shift compensation is deterministic with a fixed seed", {
  skip_if_missing_alg_deps()
  alg2_res <- alg2_rob_meanNCov(build_baseline_matrix(), print2screen = FALSE)
  source_list <- build_small_tasks()
  task_list <- build_target_tasks()

  set.seed(2024)
  first_run <- alg3_shiftComp(
    source_list = source_list,
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    print2screen = FALSE,
    save2file = FALSE,
    maximumLag = 2,
    ImpFeats = NA,
    ADM = FALSE,
    datatyp = "FC",
    useAbsCor = TRUE,
    medianMediansBL = FALSE,
    CoreClassifier = "LinSVM"
  )

  set.seed(2024)
  second_run <- alg3_shiftComp(
    source_list = source_list,
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    print2screen = FALSE,
    save2file = FALSE,
    maximumLag = 2,
    ImpFeats = NA,
    ADM = FALSE,
    datatyp = "FC",
    useAbsCor = TRUE,
    medianMediansBL = FALSE,
    CoreClassifier = "LinSVM"
  )

  expect_equal(first_run, second_run)
  expect_length(first_run, length(task_list))
  expect_true(all(is.finite(first_run)))
})

test_that("bias update aligns dimensions and reports NA inputs", {
  skip_if_missing_alg_deps()
  alg2_res <- alg2_rob_meanNCov(build_baseline_matrix(), print2screen = FALSE)
  source_list <- build_small_tasks()
  task_list <- build_target_tasks()

  set.seed(123)
  alg3_res <- alg3_shiftComp(
    source_list = source_list,
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    print2screen = FALSE,
    save2file = FALSE,
    maximumLag = 2,
    ImpFeats = NA,
    ADM = FALSE,
    datatyp = "FC",
    useAbsCor = TRUE,
    medianMediansBL = FALSE,
    CoreClassifier = "LinSVM"
  )

  bias_updates <- alg4_BiasUpdate(
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    alg3_result = alg3_res,
    goodColumns = "",
    save2file = FALSE,
    Marg = 0.5,
    alg4MinFx = "gd",
    ADM = FALSE,
    useMedian = TRUE,
    ZnormMappingBL = FALSE,
    datatyp = "FC",
    RCSmodeBL = FALSE,
    RCSfreqSet = c(0, 0),
    CoreClassifier = "LinSVM"
  )

  expect_s3_class(bias_updates, "data.frame")
  expect_equal(nrow(bias_updates), length(task_list))
  expect_false(any(is.na(bias_updates)))

  na_task <- list(task_na = data.frame(x = c(1, NA, 2)))
  expect_error(
    alg4_BiasUpdate(
      task_list = na_task,
      alg1_result = list(),
      alg2_result = alg2_res,
      alg3_result = alg3_res,
      goodColumns = "",
      save2file = FALSE,
      Marg = 0.5,
      alg4MinFx = "gd",
      ADM = FALSE,
      useMedian = TRUE,
      ZnormMappingBL = FALSE,
      datatyp = "FC",
      RCSmodeBL = FALSE,
      RCSfreqSet = c(0, 0),
      CoreClassifier = "LinSVM"
    ),
    "missing"
  )
})

test_that("normal vector update preserves dimensionality", {
  skip_if_missing_alg_deps()
  alg2_res <- alg2_rob_meanNCov(build_baseline_matrix(), print2screen = FALSE)
  source_list <- build_small_tasks()
  task_list <- build_target_tasks()

  set.seed(321)
  alg3_res <- alg3_shiftComp(
    source_list = source_list,
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    print2screen = FALSE,
    save2file = FALSE,
    maximumLag = 2,
    ImpFeats = NA,
    ADM = FALSE,
    datatyp = "FC",
    useAbsCor = TRUE,
    medianMediansBL = FALSE,
    CoreClassifier = "LinSVM"
  )

  alg4_res <- alg4_BiasUpdate(
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    alg3_result = alg3_res,
    goodColumns = "",
    save2file = FALSE,
    Marg = 0.5,
    alg4MinFx = "gd",
    ADM = FALSE,
    useMedian = TRUE,
    ZnormMappingBL = FALSE,
    datatyp = "FC",
    RCSmodeBL = FALSE,
    RCSfreqSet = c(0, 0),
    CoreClassifier = "LinSVM"
  )

  alg6_res <- alg6_NormalVectorUpdate(
    task_list = task_list,
    alg1_result = list(),
    alg2_result = alg2_res,
    alg3_result = alg3_res,
    alg4_result = alg4_res,
    X_feat_cols = "",
    Marg = 0.5,
    save2file = FALSE,
    ADM = FALSE,
    datatyp = "FC",
    RCSmodeBL = FALSE,
    RCSfreqSet = c(0, 0),
    CoreClassifier = "LinSVM"
  )

  expect_equal(nrow(alg6_res$alg6_w_new), length(task_list))
  expect_true(all(is.finite(alg6_res$alg6_w_new)))
  expect_length(alg6_res$alg6_slope, length(task_list))
})

test_that("profiling highlights hotspot routines", {
  skip_if_missing_alg_deps()
  tmp <- tempfile()
  Rprof(tmp, interval = 0.001)
  replicate(50, KBand_fx(1:3, c(1, 2, 3)))
  Rprof(NULL)

  prof <- summaryRprof(tmp)
  expect_true(file.exists(tmp))
  expect_true("KBand_fx" %in% rownames(prof$by.self))
})
