test_that("baseline classifier returns expected structure", {
  skip_if_not_installed("e1071")
  skip_if_not_installed("caret")

  train_data <- data.frame(feature = c(0, 1))
  train_labels <- factor(c("Neg", "Pos"), levels = c("Neg", "Pos"))

  train_inputs <- list(sample = train_data)
  label_inputs <- list(sample = train_labels)

  results <- alg1_baselineClass(
    TrainXls = train_inputs,
    TrainYls = label_inputs,
    TestXls = train_inputs,
    TestYls = label_inputs,
    K_forCrossV = 2,
    svmGamma = 0.1,
    svmCost = 1,
    prnt2scr = FALSE,
    X_cols2Keep = 1
  )

  expect_equal(dim(results$baselineSVM), c(1, 2))
  expect_snapshot({
    list(
      dims = dim(results$baselineSVM),
      colnames = colnames(results$baselineSVM),
      rownames = rownames(results$baselineSVM)
    )
  })
})
