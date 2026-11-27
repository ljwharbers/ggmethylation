# summary() / print() for methylation_data and multi_methylation_data.

test_that("summary works on an empty methylation_data object", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", NULL)
  expect_silent(s <- summary(md))
  expect_s3_class(s, "summary.methylation_data")
  expect_equal(s$n_reads, 0L)
  expect_equal(s$n_sites, 0L)
  expect_equal(nrow(s$strand), 0L)
  expect_named(s$strand, c("strand", "n_reads"))
  expect_true(is.na(s$read_length$median))
  out = utils::capture.output(print(s))
  expect_true(any(grepl("Reads:", out)))
  expect_true(any(grepl("no reads", out)))
  expect_false(any(grepl("min=NA", out, fixed = TRUE)))
})

test_that("summary works on an empty grouped object", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", "HP")
  expect_silent(s <- summary(md))
  expect_equal(s$n_reads, 0L)
  expect_silent(utils::capture.output(print(s)))
})

test_that("summary works when every group is NA", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", "HP")
  md$reads = data.frame(
    read_name = "r1", start = 1L, end = 10L, bam_pos = 1L, strand = "+",
    is_supplementary = FALSE, sa_chrom = NA_character_, sa_pos = NA_integer_,
    clip_side = NA_character_, group = NA_character_, stringsAsFactors = FALSE
  )
  md$sites = data.frame(position = 5L, mod_prob = 0.5, read_name = "r1",
                        mod_code = "m", group = NA_character_, stringsAsFactors = FALSE)
  s = summary(md)
  expect_equal(s$n_reads, 1L)
  expect_equal(s$overall_mean_mod_prob, 0.5)
  expect_equal(NROW(s$groups), 0L)
  # No bam_end column: the length falls back to end - start + 1.
  expect_equal(s$read_length$median, 10L)
  expect_silent(utils::capture.output(print(md)))
})

test_that("summary prints nothing itself; printing the result shows it", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000",
                                         group_tag = "HP", seed = 1))
  expect_silent(s <- summary(md))
  expect_visible(summary(md))
  expect_output(print(s), "^methylation_data summary")
  expect_output(print(s), "Aligned length: +min=[0-9]+, median=[0-9]+, max=[0-9]+ bp")
  expect_output(print(s), "Group breakdown \\(HP\\)")
  utils::capture.output(out <- withVisible(print(s)))
  expect_false(out$visible)
  expect_identical(out$value, s)
})

test_that("aligned length uses the full alignment, not the span clipped to the region", {
  # Reads run past this 201-bp window, so clipped spans are at most 201 bp.
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:3000-3200", seed = 1))
  full = md$reads$bam_end - md$reads$bam_pos + 1L
  expect_gt(max(full), 201L)
  s = summary(md)
  expect_equal(s$read_length$median, as.integer(stats::median(full)))
  expect_equal(s$read_length$max, max(full))
  expect_output(print(md), sprintf("Median aligned length: %d bp",
                                   as.integer(stats::median(full))))
})

test_that("summary of a multi_methylation_data object is a list of sample summaries", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000",
                                         group_tag = "HP", seed = 1))
  empty = empty_methylation_data(md$region, "m", NULL)
  merged = merge_methylation(A = md, B = empty)
  expect_silent(s <- summary(merged))
  expect_s3_class(s, "summary.multi_methylation_data")
  expect_named(s, c("A", "B"))
  expect_s3_class(s$A, "summary.methylation_data")
  expect_equal(s$A$n_reads, nrow(md$reads))
  expect_equal(s$B$n_reads, 0L)
  out = utils::capture.output(print(s))
  expect_true(any(grepl("=== Sample: A ===", out, fixed = TRUE)))
  expect_true(any(grepl("=== Sample: B ===", out, fixed = TRUE)))
  expect_equal(sum(grepl("^methylation_data summary", out)), 2L)
  expect_true(any(grepl("no reads", out)))
})
