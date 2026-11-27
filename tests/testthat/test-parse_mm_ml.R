# Tests for parse_mm_ml() and .parse_mods() in R/parse_mm_ml.R

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

test_that("parse_mm_ml returns empty df (with a warning) on empty ml_tag", {
  # One MM delta but no ML value: the tags disagree, so no calls are made.
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("ACGACG", "C+m,0;", integer(0), "m", "+", "6M", 1L),
    "ML"
  )
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

test_that("parse_mm_ml returns empty (with a warning) when deltas run past the sequence", {
  # seq = "ACGT": 1 C at pos 2
  # "C+m,1;" -> delta=1 -> current = 0+1+1 = 2 -> 2nd C -> out of bounds (only 1 C)
  # The MM tag does not fit SEQ (e.g. a hard-clipped record), so it is rejected.
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("ACGT", "C+m,1;", as.integer(200), "m", "+", "4M", 1L),
    "more"
  )
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

# --- MM strand field ('+' / '-') ---
# SAM spec (SAMtags, MM): the strand field only says which strand carries the
# modification. Bases are always counted on the original top-strand read
# sequence, 5' -> 3'. Only the alignment strand (FLAG 0x10) decides whether
# that sequence is SEQ as stored or its reverse complement.

test_that("'-' MM entry on a forward read counts G's in SEQ from the left", {
  # SEQ "ACGTG": G's at 3 and 5. "G-m,0" -> first G -> query 3 -> ref 102.
  result = ggmethylation:::parse_mm_ml("ACGTG", "G-m,0;", 255L, "m", "+", "5M", 100L)
  expect_equal(result$sites$position, 102L)
  expect_equal(result$sites$mod_prob, 1)
})

test_that("'-' MM entry on a reverse read counts G's in the original read", {
  # Original read (5'->3') "TGCAG": G's at 2 and 5. "G-m,1" -> 2nd G = read pos 5.
  # Stored SEQ = revcomp("TGCAG") = "CTGCA"; read pos k <-> SEQ pos 6 - k,
  # so read pos 5 -> SEQ pos 1 ('C') -> ref 100 with 5M at POS 100.
  result = ggmethylation:::parse_mm_ml("CTGCA", "G-m,1;", 128L, "m", "-", "5M", 100L)
  expect_equal(result$sites$position, 100L)
  expect_equal(result$sites$mod_prob, 128 / 255)
})

# --- Combined modification codes ---

test_that("combined codes (C+mh) read interleaved ML values", {
  # SEQ "ACGACG": C's at 2 and 5. "C+mh,0,0" lists both C's for m and h.
  # ML order: C1-m, C1-h, C2-m, C2-h = 200, 10, 100, 20.
  ml = c(200L, 10L, 100L, 20L)
  m = ggmethylation:::parse_mm_ml("ACGACG", "C+mh,0,0;", ml, "m", "+", "6M", 1L)
  expect_equal(m$sites$position, c(2L, 5L))
  expect_equal(m$sites$mod_prob, c(200, 100) / 255)
  h = ggmethylation:::parse_mm_ml("ACGACG", "C+mh,0,0;", ml, "h", "+", "6M", 1L)
  expect_equal(h$sites$position, c(2L, 5L))
  expect_equal(h$sites$mod_prob, c(10, 20) / 255)
})

test_that("a combined-code entry consumes one ML value per code per base", {
  # "C+mh,0" uses ML 200 (m) and 10 (h); "A+a,0" then uses 77.
  # SEQ "ACGAAA": first A at 1.
  result = ggmethylation:::parse_mm_ml("ACGAAA", "C+mh,0;A+a,0;", c(200L, 10L, 77L),
                                       "a", "+", "6M", 1L)
  expect_equal(result$sites$position, 1L)
  expect_equal(result$sites$mod_prob, 77 / 255)
})

test_that("a numeric ChEBI code is a single code", {
  # "C+76792,0" is one code (N4-methylcytosine), one ML value (99), then
  # "C+m,0" uses 150. SEQ "ACG": C at 2.
  chebi = ggmethylation:::parse_mm_ml("ACG", "C+76792,0;C+m,0;", c(99L, 150L),
                                      "76792", "+", "3M", 1L)
  expect_equal(chebi$sites$position, 2L)
  expect_equal(chebi$sites$mod_prob, 99 / 255)
  m = ggmethylation:::parse_mm_ml("ACG", "C+76792,0;C+m,0;", c(99L, 150L),
                                  "m", "+", "3M", 1L)
  expect_equal(m$sites$mod_prob, 150 / 255)
})

test_that("two MM entries with the same code are both returned", {
  # SEQ "ACGT": "C+m,0" -> C at 2 (ML 200); "G-m,0" -> G at 3 (ML 100).
  result = ggmethylation:::parse_mm_ml("ACGT", "C+m,0;G-m,0;", c(200L, 100L),
                                       "m", "+", "4M", 1L)
  ord = order(result$sites$position)
  expect_equal(result$sites$position[ord], c(2L, 3L))
  expect_equal(result$sites$mod_prob[ord], c(200, 100) / 255)
})

# --- Malformed tags ---

test_that("ML shorter than MM is rejected with a warning", {
  # Two deltas but one ML value.
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("ACGACG", "C+m,0,0;", 200L, "m", "+", "6M", 1L),
    "ML"
  )
  expect_equal(nrow(result$sites), 0L)
  expect_false(anyNA(result$sites$mod_prob))
})

test_that("ML longer than MM is rejected with a warning", {
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("ACG", "C+m,0;", c(200L, 100L), "m", "+", "3M", 1L),
    "ML"
  )
  expect_equal(nrow(result$sites), 0L)
})

test_that("a non-numeric MM delta is rejected with a warning, not an error", {
  expect_warning(
    r1 <- ggmethylation:::parse_mm_ml("ACGACGACG", "C+m,,1", c(1L, 2L), "m", "+", "9M", 1L),
    "delta"
  )
  expect_equal(nrow(r1$sites), 0L)
  expect_warning(
    r2 <- ggmethylation:::parse_mm_ml("ACGACGACG", "C+m,0,,1;", c(1L, 2L, 3L),
                                      "m", "+", "9M", 1L),
    "delta"
  )
  expect_equal(nrow(r2$sites), 0L)
})

test_that("an unparseable MM entry is rejected with a warning", {
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("ACG", "C*m,0;", 200L, "m", "+", "3M", 1L),
    "MM"
  )
  expect_equal(nrow(result$sites), 0L)
})

test_that("canonical base N matches any base, on either strand", {
  # Forward: "N+n,1" -> skip one base -> query 2 -> ref 2.
  fwd = ggmethylation:::parse_mm_ml("ACGT", "N+n,1;", 200L, "n", "+", "4M", 1L)
  expect_equal(fwd$sites$position, 2L)
  # Reverse: the original read starts at the end of SEQ, so "N+n,0" is
  # SEQ pos 4 -> ref 4.
  rev = ggmethylation:::parse_mm_ml("ACGT", "N+n,0;", 200L, "n", "-", "4M", 1L)
  expect_equal(rev$sites$position, 4L)
  expect_equal(rev$sites$mod_prob, 200 / 255)
})

# --- MN tag (hard clipping) ---

test_that("calls are dropped when MN differs from the SEQ length", {
  # Hard-clipped record: SEQ is 10 bases but MM was written for 15.
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("ACGTACGTAC", "C+m,0;", 200L, "m", "+",
                                          "5H10M", 1L, mn = 15L),
    "MN"
  )
  expect_equal(nrow(result$sites), 0L)
  # MN equal to the SEQ length is fine: first C at query 2 -> ref 2.
  ok = ggmethylation:::parse_mm_ml("ACGTACGTAC", "C+m,0;", 200L, "m", "+",
                                   "5H10M", 1L, mn = 10L)
  expect_equal(ok$sites$position, 2L)
})

test_that("an empty SEQ ('*') with MM calls is rejected, not crashed on", {
  expect_warning(
    result <- ggmethylation:::parse_mm_ml("", "C+m,0;", 200L, "m", "+", "4M", 1L),
    "SEQ"
  )
  expect_equal(nrow(result$sites), 0L)
})

test_that("an MM entry without deltas or ML values is valid and empty", {
  expect_silent(
    result <- ggmethylation:::parse_mm_ml("ACG", "C+m?;", integer(0), "m", "+", "3M", 1L)
  )
  expect_equal(nrow(result$sites), 0L)
})

# --- One-pass extraction of several codes ---

test_that(".parse_mods extracts every requested code in one pass", {
  # SEQ "ACGACG" (C at 2, 5), "C+m,0,0;C+h,1;": m -> 2 (200), 5 (100); h -> 5 (30).
  res = ggmethylation:::.parse_mods("ACGACG", "C+m,0,0;C+h,1;", c(200L, 100L, 30L),
                                    c("m", "h"), "+", "6M", 1L)
  expect_equal(res$sites$position, c(2L, 5L, 5L))
  expect_equal(res$sites$mod_code, c("m", "m", "h"))
  expect_equal(res$sites$mod_prob, c(200, 100, 30) / 255)
})
