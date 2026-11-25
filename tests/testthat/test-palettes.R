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

test_that(".resolve_group_colours keeps named colours and extends for others", {
  res = ggmethylation:::.resolve_group_colours
  ext = ggmethylation:::.GROUP_PALETTE_EXTENSION
  def = ggmethylation:::.GROUP_PALETTE_DEFAULT
  expect_equal(res(def, c("2", "1", NA)), def[c("1", "2")])
  # HP with a third haplotype: 1/2 keep the defaults, 3 gets the next colour
  expect_equal(res(def, c("1", "2", "3")),
               c("1" = "#0072B2", "2" = "#E69F00", "3" = ext[3L]))
  # SNV groups: none named, generated in level order, silently for the default
  expect_silent(out <- res(def, c("REF", "ALT")))
  expect_equal(out, c(ALT = ext[1L], REF = ext[2L]))
  # a user palette that misses a group says so
  expect_message(res(c(A = "red"), c("A", "B")), "does not name group")
  # unnamed: by position over the sorted levels
  expect_equal(res(c("red", "blue", "green"), c("b", "a")), c(a = "red", b = "blue"))
  expect_error(res(c("red"), c("a", "b")), "1 unnamed colour")
  # NULL: ggplot2 hue palette
  expect_equal(unname(res(NULL, c("a", "b"))), scales::hue_pal()(2))
  # no groups at all
  expect_length(res(def, c(NA, NA)), 0L)
})
