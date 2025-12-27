test_that("deterministic helper utilities remain stable", {
  helper_values <- list(
    mode_numeric = Mode(c(1, 2, 2, 3)),
    mode_character = Mode(c("a", "b", "a", "c")),
    indicator_true = IndicatorFX(TRUE),
    indicator_false = IndicatorFX(FALSE),
    odd_flags = is.odd(1:4),
    even_flags = is.even(1:4)
  )

  expect_snapshot({
    print(helper_values)
  })
})
