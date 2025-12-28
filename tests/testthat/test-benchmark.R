test_that("benchmark helper reports elapsed timings", {
  skip_if_not_installed("e1071")
  skip_if_not_installed("caret")
  skip_if_not_installed("GSE")
  skip_if_not_installed("gradDescent")

  res <- rtl_benchmark(
    sizes = c(tiny = 50),
    features = 2,
    repetitions = 1,
    margin = 0.2,
    seed = 42
  )

  expect_s3_class(res, "data.table")
  expect_true(all(res$mode %in% c("sequential", "parallel")))
  expect_true(all(res$elapsed > 0))
  expect_equal(unique(res$dataset), "tiny")
})
