# Read-level filters and downsampling of read_methylation() on the fixture.

read_all = function(...) {
  suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000", ...))
}

test_that("strand_filter keeps only the requested strand", {
  both = read_all(max_reads = Inf)
  plus = read_all(max_reads = Inf, strand_filter = "+")
  minus = read_all(max_reads = Inf, strand_filter = "-")

  expect_true(all(c("+", "-") %in% both$reads$strand))
  expect_true(all(plus$reads$strand == "+"))
  expect_true(all(minus$reads$strand == "-"))
  expect_equal(nrow(plus$reads) + nrow(minus$reads), nrow(both$reads))
  # Calls follow their reads
  expect_true(all(plus$sites$read_name %in% plus$reads$read_name))
  expect_error(read_all(strand_filter = "x"), "should be one of")
})

test_that("per_group_downsample caps every group and drops ungrouped reads", {
  full = read_all(group_tag = "HP", max_reads = Inf)
  sizes = table(full$reads$group)
  expect_true(all(sizes > 5L))
  expect_true(anyNA(full$reads$group))

  per_group = read_all(group_tag = "HP", max_reads = 5, per_group_downsample = TRUE,
                       seed = 1)
  expect_false(anyNA(per_group$reads$group))
  expect_equal(as.integer(table(per_group$reads$group)), rep(5L, length(sizes)))
  expect_setequal(unique(per_group$sites$group), names(sizes))

  # Without it the cap is global, over all reads including ungrouped ones.
  global = read_all(group_tag = "HP", max_reads = 5, seed = 1)
  expect_equal(nrow(global$reads), 5L)
})

test_that("per_group_downsample is ignored for ungrouped data", {
  md = read_all(max_reads = 7, per_group_downsample = TRUE, seed = 1)
  expect_equal(nrow(md$reads), 7L)
})
