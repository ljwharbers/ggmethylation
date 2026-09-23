test_that("theme_ggmethylation returns a ggplot2 theme", {
  th = ggmethylation::theme_ggmethylation()
  expect_s3_class(th, "theme")
})

test_that("palette constants have expected shape", {
  pg = ggmethylation:::.PROB_GRADIENT
  expect_true(all(c("low", "high") %in% names(pg)))
  gp = ggmethylation:::.GROUP_PALETTE_DEFAULT
  expect_true(all(c("1", "2") %in% names(gp)))
  dv = ggmethylation:::.DELTA_DIVERGING
  expect_true(all(c("neg", "pos", "zero") %in% names(dv)))
})

test_that("ambiguous call colour is distinguishable from the unmethylated end", {
  # Regression guard: the ambiguous colour used to be "grey75" (191,191,191),
  # two units per channel away from the default colour_low "#BDBDBD"
  # (189,189,189) - i.e. invisible as a separate category.
  amb = ggmethylation:::.CALL_AMBIGUOUS_DEFAULT
  expect_type(amb, "character")
  expect_length(amb, 1L)

  dist = sum(abs(grDevices::col2rgb(amb) - grDevices::col2rgb("#BDBDBD")))
  expect_gt(dist, 60)
})
