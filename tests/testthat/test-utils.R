# Tests for parse_region, ref_to_seq, and complement_base in R/utils.R

# --- parse_region ---

test_that("parse_region parses standard region", {
  result <- ggmethylation:::parse_region("chr1:1000-2000")
  expect_equal(result$chrom, "chr1")
  expect_equal(result$start, 1000L)
  expect_equal(result$end, 2000L)
})

test_that("parse_region strips commas", {
  result <- ggmethylation:::parse_region("chr1:1,000-2,000")
  expect_equal(result$chrom, "chr1")
  expect_equal(result$start, 1000L)
  expect_equal(result$end, 2000L)
})

test_that("parse_region rejects missing dash", {
  expect_error(ggmethylation:::parse_region("chr1:1000"), "Invalid region format")
})

test_that("parse_region rejects missing colon", {
  expect_error(ggmethylation:::parse_region("chr11000-2000"), "Invalid region format")
})

test_that("parse_region rejects empty string", {
  expect_error(ggmethylation:::parse_region(""))
})

test_that("parse_region accepts dotted chromosome names", {
  result <- ggmethylation:::parse_region("chr1.1:100-200")
  expect_equal(result$chrom, "chr1.1")
  expect_equal(result$start, 100L)
  expect_equal(result$end, 200L)
})

# --- ref_to_seq ---
# ref_to_seq(cigar, ref_start, target_ref_pos) -> single integer query position or NA

test_that("ref_to_seq M-only alignment maps correctly", {
  # 10M from ref 1000; ref pos 1000 -> query 1, ref pos 1004 -> query 5
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1000L), 1L)
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1004L), 5L)
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1009L), 10L)
})

test_that("ref_to_seq returns NA for position outside alignment", {
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 999L), NA_integer_)
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1010L), NA_integer_)
})

test_that("ref_to_seq returns NA for position in deletion", {
  # 5M3D5M from ref 100; deletion covers ref 105-107
  expect_equal(ggmethylation:::ref_to_seq("5M3D5M", 100L, 105L), NA_integer_)
  expect_equal(ggmethylation:::ref_to_seq("5M3D5M", 100L, 107L), NA_integer_)
})

test_that("ref_to_seq maps correctly after deletion", {
  # 5M3D5M from ref 100; ref 108 is first base of second 5M -> query pos 6
  expect_equal(ggmethylation:::ref_to_seq("5M3D5M", 100L, 108L), 6L)
})

test_that("ref_to_seq soft clip shifts effective query positions", {
  # 4S6M from ref 100; ref 100 -> query 5, ref 105 -> query 10
  expect_equal(ggmethylation:::ref_to_seq("4S6M", 100L, 100L), 5L)
  expect_equal(ggmethylation:::ref_to_seq("4S6M", 100L, 105L), 10L)
  # Soft-clipped region has no reference position -> query 1-4 are not reachable by ref coord
  expect_equal(ggmethylation:::ref_to_seq("4S6M", 100L, 99L), NA_integer_)
})

# --- complement_base ---

test_that("complement_base returns correct complements", {
  expect_equal(ggmethylation:::complement_base("A"), "T")
  expect_equal(ggmethylation:::complement_base("T"), "A")
  expect_equal(ggmethylation:::complement_base("C"), "G")
  expect_equal(ggmethylation:::complement_base("G"), "C")
})

test_that("complement_base errors on unknown base", {
  expect_error(ggmethylation:::complement_base("N"), "Unknown base")
})
