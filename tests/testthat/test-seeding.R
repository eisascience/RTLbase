test_that("scoped_seed restores RNG state after use", {
  set.seed(999)
  baseline <- runif(1)

  set.seed(999)
  scoped_seed(123, runif(5))

  after <- runif(1)
  expect_equal(baseline, after)
})

test_that("scoped seeds enable snapshot-friendly artifacts", {
  seeded_df <- scoped_seed(1, data.frame(x = runif(5), y = runif(5)))
  plot_obj <- ggplot2::ggplot(seeded_df, ggplot2::aes(x, y)) + ggplot2::geom_point()

  expect_snapshot(seeded_df)
  expect_snapshot(ggplot2::ggplot_build(plot_obj)$data[[1]][c("x", "y")])
})

