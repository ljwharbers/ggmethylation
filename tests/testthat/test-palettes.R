test_that("theme_ggmethylation returns a ggplot2 theme", {
  th <- ggmethylation::theme_ggmethylation()
  expect_s3_class(th, "theme")
})

test_that("palette constants have expected shape", {
  expect_true(is.character(ggmethylation:::.OKABE_ITO))
  expect_gte(length(ggmethylation:::.OKABE_ITO), 8L)
  gp <- ggmethylation:::.GROUP_PALETTE_DEFAULT
  expect_true(all(c("1", "2") %in% names(gp)))
  dv <- ggmethylation:::.DELTA_DIVERGING
  expect_true(all(c("neg", "pos", "zero") %in% names(dv)))
})
