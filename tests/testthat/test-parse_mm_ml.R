# Tests for seq_to_ref and parse_mm_ml in R/parse_mm_ml.R

# --- seq_to_ref ---
# seq_to_ref(cigar, pos, query_positions) -> ref positions for each query position

test_that("seq_to_ref M-only alignment maps query positions to ref", {
  # 10M from ref 1000: query pos 1 -> ref 1000, pos 5 -> ref 1004, pos 10 -> ref 1009
  result <- ggmethylation:::seq_to_ref("10M", 1000L, c(1L, 5L, 10L))
  expect_equal(result, c(1000L, 1004L, 1009L))
})

test_that("seq_to_ref soft clip returns NA for clipped query positions", {
  # 4S6M from ref 100: query pos 1-4 are soft-clipped -> NA
  # query pos 5 -> ref 100, pos 10 -> ref 105
  result <- ggmethylation:::seq_to_ref("4S6M", 100L, c(1L, 4L, 5L, 10L))
  expect_equal(result, c(NA_integer_, NA_integer_, 100L, 105L))
})

test_that("seq_to_ref deletion skips reference positions", {
  # 5M3D5M from ref 100:
  #   query 1-5 -> ref 100-104
  #   deletion skips ref 105-107
  #   query 6-10 -> ref 108-112
  result <- ggmethylation:::seq_to_ref("5M3D5M", 100L, c(5L, 6L))
  expect_equal(result, c(104L, 108L))
})

test_that("seq_to_ref insertion positions return NA", {
  # 3M2I4M from ref 1:
  #   query 1-3 -> ref 1-3
  #   query 4-5 are insertion -> NA
  #   query 6-9 -> ref 4-7
  result <- ggmethylation:::seq_to_ref("3M2I4M", 1L, c(3L, 4L, 5L, 6L))
  expect_equal(result, c(3L, NA_integer_, NA_integer_, 4L))
})

test_that("seq_to_ref returns same length as query_positions", {
  result <- ggmethylation:::seq_to_ref("10M", 1L, c(1L, 3L, 7L, 10L))
  expect_equal(length(result), 4L)
})

# --- parse_mm_ml ---

test_that("parse_mm_ml forward strand simple CpG", {
  # seq = "ACGACGACG", 3 C's at positions 2, 5, 8
  # MM: "C+m,0,1;" -> delta=0 -> C index 1 (pos 2), delta=1 -> C index 3 (pos 8)
  # (delta 0 = next C is this one; delta 1 = skip 1 C, so C index 1+1+1=3)
  # With pos=1000, cigar="9M": ref pos = query pos + 999
  # C at query 2 -> ref 1001; C at query 8 -> ref 1007
  seq    <- "ACGACGACG"
  mm_tag <- "C+m,0,1;"
  ml_tag <- as.integer(c(230, 128))
  result <- ggmethylation:::parse_mm_ml(seq, mm_tag, ml_tag, "m", "+", "9M", 1000L)
  expect_equal(result$position, c(1001L, 1007L))
  expect_equal(result$mod_prob, c(230 / 255, 128 / 255))
})

test_that("parse_mm_ml returns empty df when mod_code absent", {
  result <- ggmethylation:::parse_mm_ml("ACGACG", "C+m,0;", as.integer(c(200)), "h", "+", "6M", 1L)
  expect_equal(nrow(result), 0L)
  expect_true("position" %in% names(result))
  expect_true("mod_prob" %in% names(result))
})

test_that("parse_mm_ml returns empty df on NULL mm_tag", {
  result <- ggmethylation:::parse_mm_ml("ACGACG", NULL, as.integer(c(200)), "m", "+", "6M", 1L)
  expect_equal(nrow(result), 0L)
})

test_that("parse_mm_ml returns empty df on NA mm_tag", {
  result <- ggmethylation:::parse_mm_ml("ACGACG", NA_character_, as.integer(c(200)), "m", "+", "6M", 1L)
  expect_equal(nrow(result), 0L)
})

test_that("parse_mm_ml returns empty df on empty ml_tag", {
  result <- ggmethylation:::parse_mm_ml("ACGACG", "C+m,0;", integer(0), "m", "+", "6M", 1L)
  expect_equal(nrow(result), 0L)
})

test_that("parse_mm_ml mod_prob is ml_value / 255", {
  # seq "ACG": C at pos 2; "C+m,0;" -> delta=0 -> first C
  result <- ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(255), "m", "+", "3M", 1L)
  expect_equal(result$mod_prob, 1.0)
  result2 <- ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(0), "m", "+", "3M", 1L)
  expect_equal(result2$mod_prob, 0.0)
})

test_that("parse_mm_ml handles trailing semicolon vs no semicolon identically", {
  r1 <- ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(200), "m", "+", "3M", 1L)
  r2 <- ggmethylation:::parse_mm_ml("ACG", "C+m,0", as.integer(200), "m", "+", "3M", 1L)
  expect_equal(r1$position, r2$position)
  expect_equal(r1$mod_prob, r2$mod_prob)
})

test_that("parse_mm_ml returns empty when all deltas out of bounds", {
  # seq = "ACGT": 1 C at pos 2
  # "C+m,1;" -> delta=1 -> current = 0+1+1 = 2 -> 2nd C -> out of bounds (only 1 C)
  result <- ggmethylation:::parse_mm_ml("ACGT", "C+m,1;", as.integer(200), "m", "+", "4M", 1L)
  expect_equal(nrow(result), 0L)
})

test_that("parse_mm_ml result columns are position (integer) and mod_prob (numeric)", {
  result <- ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(200), "m", "+", "3M", 1L)
  expect_true(is.integer(result$position))
  expect_true(is.numeric(result$mod_prob))
})
