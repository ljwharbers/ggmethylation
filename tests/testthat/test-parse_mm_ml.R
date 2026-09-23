# Tests for seq_to_ref and parse_mm_ml in R/parse_mm_ml.R

# --- seq_to_ref ---
# seq_to_ref(cigar, pos, query_positions) -> ref positions for each query position

test_that("seq_to_ref M-only alignment maps query positions to ref", {
  # 10M from ref 1000: query pos 1 -> ref 1000, pos 5 -> ref 1004, pos 10 -> ref 1009
  result = ggmethylation:::seq_to_ref("10M", 1000L, c(1L, 5L, 10L))
  expect_equal(result, c(1000L, 1004L, 1009L))
})

test_that("seq_to_ref soft clip returns NA for clipped query positions", {
  # 4S6M from ref 100: query pos 1-4 are soft-clipped -> NA
  # query pos 5 -> ref 100, pos 10 -> ref 105
  result = ggmethylation:::seq_to_ref("4S6M", 100L, c(1L, 4L, 5L, 10L))
  expect_equal(result, c(NA_integer_, NA_integer_, 100L, 105L))
})

test_that("seq_to_ref deletion skips reference positions", {
  # 5M3D5M from ref 100:
  #   query 1-5 -> ref 100-104
  #   deletion skips ref 105-107
  #   query 6-10 -> ref 108-112
  result = ggmethylation:::seq_to_ref("5M3D5M", 100L, c(5L, 6L))
  expect_equal(result, c(104L, 108L))
})

test_that("seq_to_ref insertion positions return NA", {
  # 3M2I4M from ref 1:
  #   query 1-3 -> ref 1-3
  #   query 4-5 are insertion -> NA
  #   query 6-9 -> ref 4-7
  result = ggmethylation:::seq_to_ref("3M2I4M", 1L, c(3L, 4L, 5L, 6L))
  expect_equal(result, c(3L, NA_integer_, NA_integer_, 4L))
})

test_that("seq_to_ref returns same length as query_positions", {
  result = ggmethylation:::seq_to_ref("10M", 1L, c(1L, 3L, 7L, 10L))
  expect_equal(length(result), 4L)
})

# --- parse_mm_ml ---

test_that("parse_mm_ml forward strand simple CpG", {
  # seq = "ACGACGACG", 3 C's at positions 2, 5, 8
  # MM: "C+m,0,1;" -> delta=0 -> C index 1 (pos 2), delta=1 -> C index 3 (pos 8)
  # (delta 0 = next C is this one; delta 1 = skip 1 C, so C index 1+1+1=3)
  # With pos=1000, cigar="9M": ref pos = query pos + 999
  # C at query 2 -> ref 1001; C at query 8 -> ref 1007
  seq    = "ACGACGACG"
  mm_tag = "C+m,0,1;"
  ml_tag = as.integer(c(230, 128))
  result = ggmethylation:::parse_mm_ml(seq, mm_tag, ml_tag, "m", "+", "9M", 1000L)
  expect_equal(result$sites$position, c(1001L, 1007L))
  expect_equal(result$sites$mod_prob, c(230 / 255, 128 / 255))
})

test_that("parse_mm_ml returns a 2-element named list", {
  result = ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(200), "m", "+", "3M", 1L)
  expect_true(is.list(result))
  expect_setequal(names(result), c("sites", "insertion_sites"))
})

test_that("parse_mm_ml returns empty df when mod_code absent", {
  result = ggmethylation:::parse_mm_ml("ACGACG", "C+m,0;", as.integer(c(200)), "h", "+", "6M", 1L)
  expect_equal(nrow(result$sites), 0L)
  expect_true("position" %in% names(result$sites))
  expect_true("mod_prob" %in% names(result$sites))
})

test_that("parse_mm_ml returns empty df on NULL mm_tag", {
  result = ggmethylation:::parse_mm_ml("ACGACG", NULL, as.integer(c(200)), "m", "+", "6M", 1L)
  expect_equal(nrow(result$sites), 0L)
})

test_that("parse_mm_ml returns empty df on NA mm_tag", {
  result = ggmethylation:::parse_mm_ml("ACGACG", NA_character_, as.integer(c(200)), "m", "+", "6M", 1L)
  expect_equal(nrow(result$sites), 0L)
})

test_that("parse_mm_ml returns empty df on empty ml_tag", {
  result = ggmethylation:::parse_mm_ml("ACGACG", "C+m,0;", integer(0), "m", "+", "6M", 1L)
  expect_equal(nrow(result$sites), 0L)
})

test_that("parse_mm_ml mod_prob is ml_value / 255", {
  # seq "ACG": C at pos 2; "C+m,0;" -> delta=0 -> first C
  result = ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(255), "m", "+", "3M", 1L)
  expect_equal(result$sites$mod_prob, 1.0)
  result2 = ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(0), "m", "+", "3M", 1L)
  expect_equal(result2$sites$mod_prob, 0.0)
})

test_that("parse_mm_ml handles trailing semicolon vs no semicolon identically", {
  r1 = ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(200), "m", "+", "3M", 1L)
  r2 = ggmethylation:::parse_mm_ml("ACG", "C+m,0", as.integer(200), "m", "+", "3M", 1L)
  expect_equal(r1$sites$position, r2$sites$position)
  expect_equal(r1$sites$mod_prob, r2$sites$mod_prob)
})

test_that("parse_mm_ml returns empty when all deltas out of bounds", {
  # seq = "ACGT": 1 C at pos 2
  # "C+m,1;" -> delta=1 -> current = 0+1+1 = 2 -> 2nd C -> out of bounds (only 1 C)
  result = ggmethylation:::parse_mm_ml("ACGT", "C+m,1;", as.integer(200), "m", "+", "4M", 1L)
  expect_equal(nrow(result$sites), 0L)
})

test_that("parse_mm_ml result columns are position (integer) and mod_prob (numeric)", {
  result = ggmethylation:::parse_mm_ml("ACG", "C+m,0;", as.integer(200), "m", "+", "3M", 1L)
  expect_true(is.integer(result$sites$position))
  expect_true(is.numeric(result$sites$mod_prob))
})

test_that("parse_mm_ml routes insertion-region mods to $insertion_sites", {
  # seq: 13 bases, CIGAR 5M3I5M from pos 100
  #   query 1-5  -> ref 100-104  (M)
  #   query 6-8  -> NA, insertion at ref_anchor 105  (I)
  #   query 9-13 -> ref 105-109  (M)
  # C's at query positions 2 (M), 7 (I), 12 (M)
  # MM "C+m,0,0,0;" picks all 3; ML = 200, 100, 50
  result = ggmethylation:::parse_mm_ml(
    seq      = "ACAAAACAAAACA",
    mm_tag   = "C+m,0,0,0;",
    ml_tag   = as.integer(c(200, 100, 50)),
    mod_code = "m",
    strand   = "+",
    cigar    = "5M3I5M",
    pos      = 100L
  )
  expect_equal(result$sites$position, c(101L, 108L))
  expect_equal(result$sites$mod_prob, c(200 / 255, 50 / 255))
  expect_equal(result$insertion_sites$query_pos, 7L)
  expect_equal(result$insertion_sites$mod_prob, 100 / 255)
})

test_that("parse_mm_ml drops mods that fall in soft clips", {
  # seq: 8 bases, CIGAR 3S5M from pos 100
  #   query 1-3 -> NA (soft clip)
  #   query 4-8 -> ref 100-104
  # C's at query 2 (soft clip) and query 6 (M -> ref 102)
  # MM "C+m,0,0;" picks both; ML = 200, 100
  result = ggmethylation:::parse_mm_ml(
    seq      = "ACAAACAA",
    mm_tag   = "C+m,0,0;",
    ml_tag   = as.integer(c(200, 100)),
    mod_code = "m",
    strand   = "+",
    cigar    = "3S5M",
    pos      = 100L
  )
  expect_equal(result$sites$position, 102L)
  expect_equal(nrow(result$insertion_sites), 0L)
})

test_that("parse_mm_ml routes insertion mods to $insertion_sites on minus strand", {
  # seq: 13 bases (BAM SEQ), CIGAR 5M3I5M from pos 100
  # G at query positions 3 (M: ref 102), 7 (I: insertion), 12 (M: ref 108)
  # Strand "-", MM "C+m,0,0;" -> search_base="G", reverse_scan=TRUE
  # canonical_positions reversed = c(12, 7, 3); deltas 0,0 -> positions c(12, 7)
  # query 12 -> ref 108 (site); query 7 -> NA insertion -> insertion_site
  result = ggmethylation:::parse_mm_ml(
    seq      = "AAGAAAGAAAAGA",
    mm_tag   = "C+m,0,0;",
    ml_tag   = as.integer(c(200, 50)),
    mod_code = "m",
    strand   = "-",
    cigar    = "5M3I5M",
    pos      = 100L
  )
  expect_equal(result$sites$position, 108L)
  expect_equal(result$sites$mod_prob, 200 / 255)
  expect_equal(result$insertion_sites$query_pos, 7L)
  expect_equal(result$insertion_sites$mod_prob, 50 / 255)
})

test_that("parse_mm_ml routes insertion mods to $insertion_sites on minus strand", {
  # seq: 13 bases (BAM SEQ), CIGAR 5M3I5M from pos 100
  # G at query positions 3 (M: ref 102), 7 (I: insertion), 12 (M: ref 108)
  # Strand "-", MM "C+m,0,0;" -> search_base="G", reverse_scan=TRUE
  # canonical_positions reversed = c(12, 7, 3); deltas 0,0 -> positions c(12, 7)
  # query 12 -> ref 108 (site); query 7 -> NA insertion
  result = ggmethylation:::parse_mm_ml(
    seq      = "AAGAAAGAAAAGA",
    mm_tag   = "C+m,0,0;",
    ml_tag   = as.integer(c(200, 50)),
    mod_code = "m",
    strand   = "-",
    cigar    = "5M3I5M",
    pos      = 100L
  )
  expect_equal(result$sites$position, 108L)
  expect_equal(result$sites$mod_prob, 200 / 255)
  expect_equal(result$insertion_sites$query_pos, 7L)
  expect_equal(result$insertion_sites$mod_prob, 50 / 255)
})

# --- MM tag flag semantics ('.' vs '?' vs absent) ---

test_that("parse_mm_ml '.' flag emits unlisted CpGs as mod_prob = 0 (forward strand)", {
  # seq = "ACGACGACG": C's at query positions 2, 5, 8 -> ref 1001, 1004, 1007
  # "C+m.,0;" flag='.', delta=0 -> first C (pos 2, ref 1001) is listed/modified
  # Unlisted C's at pos 5 and 8 are implicitly unmodified -> mod_prob = 0
  result = ggmethylation:::parse_mm_ml(
    seq      = "ACGACGACG",
    mm_tag   = "C+m.,0;",
    ml_tag   = as.integer(200),
    mod_code = "m",
    strand   = "+",
    cigar    = "9M",
    pos      = 1000L
  )
  expect_equal(sort(result$sites$position), c(1001L, 1004L, 1007L))
  listed_row = result$sites[result$sites$position == 1001L, ]
  expect_equal(listed_row$mod_prob, 200 / 255)
  zero_rows = result$sites[result$sites$position %in% c(1004L, 1007L), ]
  expect_true(all(zero_rows$mod_prob == 0))
  expect_equal(nrow(result$sites), 3L)
})

test_that("parse_mm_ml '?' flag does NOT emit unlisted CpGs", {
  # Same setup but flag='?' -> only the listed C is emitted
  result = ggmethylation:::parse_mm_ml(
    seq      = "ACGACGACG",
    mm_tag   = "C+m?,0;",
    ml_tag   = as.integer(200),
    mod_code = "m",
    strand   = "+",
    cigar    = "9M",
    pos      = 1000L
  )
  expect_equal(result$sites$position, 1001L)
  expect_equal(result$sites$mod_prob, 200 / 255)
  expect_equal(nrow(result$sites), 1L)
})

test_that("parse_mm_ml no flag is treated conservatively (like '?')", {
  # "C+m,0;" has no flag -> only the listed C is emitted
  result = ggmethylation:::parse_mm_ml(
    seq      = "ACGACGACG",
    mm_tag   = "C+m,0;",
    ml_tag   = as.integer(200),
    mod_code = "m",
    strand   = "+",
    cigar    = "9M",
    pos      = 1000L
  )
  expect_equal(result$sites$position, 1001L)
  expect_equal(nrow(result$sites), 1L)
})

test_that("parse_mm_ml '.' flag emits implicit zeros on reverse strand", {
  # seq = "AGCAGCAGC" (BAM sequence for a '-' strand read)
  # G's at query positions 2, 5, 8
  # strand="-", "C+m.,0;" -> search_base="G", reverse_scan=TRUE
  # canonical_positions reversed: c(8, 5, 2)
  # delta=0 -> first in reversed order -> query 8 (listed, ref 1007)
  # unwalked -> query 5 (ref 1004), query 2 (ref 1001) -> mod_prob = 0
  result = ggmethylation:::parse_mm_ml(
    seq      = "AGCAGCAGC",
    mm_tag   = "C+m.,0;",
    ml_tag   = as.integer(200),
    mod_code = "m",
    strand   = "-",
    cigar    = "9M",
    pos      = 1000L
  )
  expect_equal(sort(result$sites$position), c(1001L, 1004L, 1007L))
  expect_equal(result$sites$mod_prob[result$sites$position == 1007L], 200 / 255)
  expect_true(all(result$sites$mod_prob[result$sites$position %in% c(1001L, 1004L)] == 0))
  expect_equal(nrow(result$sites), 3L)
})
