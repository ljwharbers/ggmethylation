test_that("decompose_cigar handles pure M alignment", {
  result <- decompose_cigar("10M", 100L)
  expect_equal(nrow(result), 1L)
  expect_equal(result$type, "M")
  expect_equal(result$ref_start, 100L)
  expect_equal(result$ref_end, 109L)
  expect_equal(result$query_start, 1L)
  expect_equal(result$query_end, 10L)
  expect_equal(result$length, 10L)
})

test_that("decompose_cigar handles deletion", {
  result <- decompose_cigar("5M3D5M", 100L)
  expect_equal(nrow(result), 3L)
  expect_equal(result$type, c("M", "D", "M"))

  del <- result[result$type == "D", ]
  expect_equal(del$ref_start, 105L)
  expect_equal(del$ref_end, 107L)
  expect_true(is.na(del$query_start))
  expect_true(is.na(del$query_end))
  expect_equal(del$length, 3L)

  # Second M block should start after the deletion
  m2 <- result[3, ]
  expect_equal(m2$ref_start, 108L)
  expect_equal(m2$ref_end, 112L)
  expect_equal(m2$query_start, 6L)
  expect_equal(m2$query_end, 10L)
})

test_that("decompose_cigar handles insertion", {
  result <- decompose_cigar("5M2I5M", 100L)
  expect_equal(nrow(result), 3L)
  expect_equal(result$type, c("M", "I", "M"))

  ins <- result[result$type == "I", ]
  expect_equal(ins$ref_start, 105L)  # insertion point
  expect_true(is.na(ins$ref_end))
  expect_equal(ins$query_start, 6L)
  expect_equal(ins$query_end, 7L)
  expect_equal(ins$length, 2L)

  # Second M: ref continues from 105 (insertion doesn't consume ref)
  m2 <- result[3, ]
  expect_equal(m2$ref_start, 105L)
  expect_equal(m2$query_start, 8L)
})

test_that("decompose_cigar handles soft clip", {
  result <- decompose_cigar("3S7M", 100L)
  expect_equal(nrow(result), 2L)
  expect_equal(result$type, c("S", "M"))

  sc <- result[result$type == "S", ]
  expect_equal(sc$length, 3L)
  expect_true(is.na(sc$ref_start))
  expect_true(is.na(sc$ref_end))
  expect_equal(sc$query_start, 1L)
  expect_equal(sc$query_end, 3L)

  # M block starts at query position 4
  m <- result[result$type == "M", ]
  expect_equal(m$ref_start, 100L)
  expect_equal(m$query_start, 4L)
})

test_that("decompose_cigar handles hard clip", {
  result <- decompose_cigar("5H10M5H", 100L)
  expect_equal(nrow(result), 3L)
  expect_equal(result$type, c("H", "M", "H"))

  h1 <- result[1, ]
  expect_true(is.na(h1$ref_start))
  expect_true(is.na(h1$query_start))
  expect_equal(h1$length, 5L)
})

test_that("decompose_cigar handles complex CIGAR", {
  result <- decompose_cigar("5S10M2I3M4D5M3S", 1000L)
  expect_equal(nrow(result), 7L)
  expect_equal(result$type, c("S", "M", "I", "M", "D", "M", "S"))

  # Verify ref coordinates chain correctly
  # S: no ref
  # M: 1000-1009 (10 bases)
  # I: ref_start=1010, no ref_end

  # M: 1010-1012 (3 bases)
  # D: 1013-1016 (4 bases)
  # M: 1017-1021 (5 bases)
  # S: no ref
  expect_equal(result$ref_start[2], 1000L)
  expect_equal(result$ref_end[2], 1009L)
  expect_equal(result$ref_start[3], 1010L)  # insertion point
  expect_equal(result$ref_start[4], 1010L)
  expect_equal(result$ref_end[4], 1012L)
  expect_equal(result$ref_start[5], 1013L)
  expect_equal(result$ref_end[5], 1016L)
  expect_equal(result$ref_start[6], 1017L)
  expect_equal(result$ref_end[6], 1021L)
})

test_that("decompose_cigar handles = and X operations", {
  result <- decompose_cigar("5=3X2M", 100L)
  expect_equal(nrow(result), 3L)
  expect_equal(result$type, c("=", "X", "M"))
  expect_equal(result$ref_start, c(100L, 105L, 108L))
  expect_equal(result$ref_end, c(104L, 107L, 109L))
  expect_equal(result$query_start, c(1L, 6L, 9L))
})

test_that("decompose_cigar handles N (skipped region)", {
  result <- decompose_cigar("5M100N5M", 100L)
  expect_equal(nrow(result), 3L)
  n_op <- result[result$type == "N", ]
  expect_equal(n_op$ref_start, 105L)
  expect_equal(n_op$ref_end, 204L)
  expect_true(is.na(n_op$query_start))
  expect_equal(n_op$length, 100L)
})

# --- split_cigar ---

test_that("split_cigar returns correct ops and lens for simple CIGAR", {
  result <- ggmethylation:::split_cigar("10M")
  expect_equal(result$ops,  "M")
  expect_equal(result$lens, 10L)
})

test_that("split_cigar handles multi-op CIGAR", {
  result <- ggmethylation:::split_cigar("5S10M2I3M")
  expect_equal(result$ops,  c("S", "M", "I", "M"))
  expect_equal(result$lens, c(5L, 10L, 2L, 3L))
})

test_that("split_cigar handles = and X ops", {
  result <- ggmethylation:::split_cigar("5=3X")
  expect_equal(result$ops,  c("=", "X"))
  expect_equal(result$lens, c(5L, 3L))
})
