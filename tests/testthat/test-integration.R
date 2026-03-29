# Integration tests using the real modBAM file.
# All tests are skipped automatically when the BAM is not present.

bam_path <- "/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam"

test_that("read_methylation returns methylation_data for real BAM (ungrouped)", {
  skip_if_not(file.exists(bam_path), "Real BAM not available")
  result <- read_methylation(bam_path, "chr21:10000000-10010000")
  expect_s3_class(result, "methylation_data")
  expect_true(is.data.frame(result$reads))
  expect_true(is.data.frame(result$sites))
  expect_true(all(c("read_name", "start", "end", "strand") %in% names(result$reads)))
  expect_true(all(c("position", "mod_prob", "read_name", "mod_code") %in% names(result$sites)))
  expect_true(all(result$sites$mod_prob >= 0 & result$sites$mod_prob <= 1, na.rm = TRUE))
})

test_that("read_methylation returns methylation_data with HP grouping", {
  skip_if_not(file.exists(bam_path), "Real BAM not available")
  result <- read_methylation(bam_path, "chr21:10000000-10010000", group_tag = "HP")
  expect_s3_class(result, "methylation_data")
  expect_true("group" %in% names(result$reads))
})

test_that("plot_methylation returns a ggplot/patchwork for real data", {
  skip_if_not(file.exists(bam_path), "Real BAM not available")
  md <- read_methylation(bam_path, "chr21:10000000-10010000")
  p <- plot_methylation(md)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})
