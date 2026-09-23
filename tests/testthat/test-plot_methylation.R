# plot_methylation() argument handling, on the synthetic fixture BAM.

fixture_md = function(...) {
  set.seed(1)
  suppressWarnings(read_methylation(
    testthat::test_path("fixtures", "synthetic.bam"), "chr1:1000-7000", ...
  ))
}

test_that("sort_by rejects columns that don't exist in data$reads", {
  # order() on a NULL column returns integer(0), which used to silently drop
  # every read (the vignette once passed sort_by = "position").
  md = fixture_md(group_tag = "HP")
  expect_error(plot_methylation(md, sort_by = "position"), "position")
  expect_error(
    plot_methylation(merge_methylation(A = md, B = md), sort_by = "position"),
    "position"
  )
  expect_s3_class(plot_methylation(md, sort_by = "mean_mod_prob"), "patchwork")
})
