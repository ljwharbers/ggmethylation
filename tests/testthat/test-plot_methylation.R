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

test_that("multi-sample plots mix grouped and ungrouped samples", {
  # Regression: combining the samples' site tables used to fail in rbind()
  # because only the grouped sample carried a `group` column.
  grouped   = fixture_md(group_tag = "HP")
  ungrouped = fixture_md()
  p = plot_methylation(merge_methylation(A = grouped, B = ungrouped))
  expect_s3_class(p, "patchwork")
  smooth = ggplot2::ggplot_build(p[[length(p)]])$data
  expect_length(smooth, 2L)  # CI ribbon + lines
  expect_gt(nrow(smooth[[2L]]), 0L)
})

test_that("multi-sample read panels use the same lane packing as single-sample", {
  # The multi-sample path used to skip clip-aware packing.
  md = fixture_md()
  single = plot_methylation(md)
  multi = plot_methylation(merge_methylation(A = md, B = md))
  lanes = function(p) sort(unique(ggplot2::ggplot_build(p)$data[[1L]]$y))
  expect_equal(lanes(multi[[1L]]), lanes(single[[1L]]))
})
