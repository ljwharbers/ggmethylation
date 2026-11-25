# summary() / print() for methylation_data.

test_that("summary works on an empty methylation_data object", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", NULL)
  out = utils::capture.output(s <- summary(md))
  expect_equal(s$n_reads, 0L)
  expect_equal(s$n_sites, 0L)
  expect_equal(nrow(s$strand), 0L)
  expect_named(s$strand, c("strand", "n_reads"))
  expect_true(any(grepl("Reads:", out)))
})

test_that("summary works on an empty grouped object", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", "HP")
  utils::capture.output(s <- summary(md))
  expect_equal(s$n_reads, 0L)
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
  utils::capture.output(s <- summary(md))
  expect_equal(s$n_reads, 1L)
  expect_equal(s$overall_mean_mod_prob, 0.5)
  expect_equal(NROW(s$groups), 0L)
  expect_silent(utils::capture.output(print(md)))
})
