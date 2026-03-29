# Tests for pack_reads in R/pack_reads.R
# Signature: pack_reads(reads, gap = 10)
# reads: data.frame with columns read_name, start, end
# Returns: integer vector of lane assignments (same length as nrow(reads))

test_that("pack_reads returns integer(0) for empty input", {
  reads <- data.frame(
    read_name = character(0),
    start = integer(0),
    end = integer(0)
  )
  result <- ggmethylation:::pack_reads(reads)
  expect_equal(result, integer(0))
})

test_that("pack_reads assigns lane 1 to single read", {
  reads <- data.frame(read_name = "r1", start = 100L, end = 200L)
  result <- ggmethylation:::pack_reads(reads)
  expect_equal(result, 1L)
})

test_that("pack_reads packs non-overlapping reads into lane 1", {
  # With gap=10, reads must have start > prev_end + 10 to fit same lane
  # r1 ends 200; r2 starts 300 -> 300 > 200+10=210 -> same lane
  # r2 ends 400; r3 starts 500 -> 500 > 400+10=410 -> same lane
  reads <- data.frame(
    read_name = c("r1", "r2", "r3"),
    start = c(100L, 300L, 500L),
    end   = c(200L, 400L, 600L)
  )
  result <- ggmethylation:::pack_reads(reads, gap = 10L)
  expect_equal(result, c(1L, 1L, 1L))
})

test_that("pack_reads assigns overlapping reads to different lanes", {
  # r1 ends 300; r2 starts 200 -> 200 > 300+10=310? No -> different lane
  reads <- data.frame(
    read_name = c("r1", "r2"),
    start = c(100L, 200L),
    end   = c(300L, 400L)
  )
  result <- ggmethylation:::pack_reads(reads, gap = 10L)
  expect_equal(result, c(1L, 2L))
})

test_that("pack_reads respects gap parameter", {
  # r1 ends 200; r2 starts 205
  # gap=10: 205 > 200+10=210? No -> different lanes
  reads <- data.frame(
    read_name = c("r1", "r2"),
    start = c(100L, 205L),
    end   = c(200L, 300L)
  )
  result_gap10 <- ggmethylation:::pack_reads(reads, gap = 10L)
  expect_equal(result_gap10, c(1L, 2L))

  # gap=0: 205 > 200+0=200? Yes -> same lane
  result_gap0 <- ggmethylation:::pack_reads(reads, gap = 0L)
  expect_equal(result_gap0, c(1L, 1L))
})

test_that("pack_reads result length matches nrow(reads)", {
  reads <- data.frame(
    read_name = paste0("r", 1:10),
    start = seq(1L, 91L, by = 10L),
    end   = seq(9L, 99L, by = 10L)
  )
  result <- ggmethylation:::pack_reads(reads)
  expect_equal(length(result), nrow(reads))
})

test_that("pack_reads returns integer vector", {
  reads <- data.frame(read_name = "r1", start = 1L, end = 10L)
  result <- ggmethylation:::pack_reads(reads)
  expect_true(is.integer(result))
})

test_that("pack_reads lane assignments are positive integers", {
  reads <- data.frame(
    read_name = c("r1", "r2", "r3"),
    start = c(1L, 5L, 9L),
    end   = c(6L, 10L, 14L)
  )
  result <- ggmethylation:::pack_reads(reads, gap = 0L)
  expect_true(all(result >= 1L))
})
