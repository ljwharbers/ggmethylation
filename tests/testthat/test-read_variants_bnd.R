# Tests for parse_bnd_alt() in R/read_variants.R
# parse_bnd_alt(alt_str) -> list(mate_chrom, mate_pos)  or NA for non-BND strings.

# --- Non-BND inputs return NA --------------------------------------------------

test_that("parse_bnd_alt returns NA for a plain SNV ALT", {
  result <- ggmethylation:::parse_bnd_alt("A")
  expect_equal(result$mate_chrom, NA_character_)
  expect_equal(result$mate_pos,   NA_integer_)
})

test_that("parse_bnd_alt returns NA for a symbolic SV ALT", {
  result <- ggmethylation:::parse_bnd_alt("<DEL>")
  expect_equal(result$mate_chrom, NA_character_)
  expect_equal(result$mate_pos,   NA_integer_)
})

test_that("parse_bnd_alt returns NA for NA input", {
  result <- ggmethylation:::parse_bnd_alt(NA_character_)
  expect_equal(result$mate_chrom, NA_character_)
  expect_equal(result$mate_pos,   NA_integer_)
})

test_that("parse_bnd_alt returns NA for empty string", {
  result <- ggmethylation:::parse_bnd_alt("")
  expect_equal(result$mate_chrom, NA_character_)
  expect_equal(result$mate_pos,   NA_integer_)
})

# --- Orientation 1: N[chr:pos[  (nucleotide, open bracket, mate, open bracket) --

test_that("parse_bnd_alt handles N[chr:pos[ orientation", {
  result <- ggmethylation:::parse_bnd_alt("A[chr2:12345[")
  expect_equal(result$mate_chrom, "chr2")
  expect_equal(result$mate_pos,   12345L)
})

test_that("parse_bnd_alt handles dot.[chr:pos[ orientation", {
  result <- ggmethylation:::parse_bnd_alt(".[chr3:99999[")
  expect_equal(result$mate_chrom, "chr3")
  expect_equal(result$mate_pos,   99999L)
})

# --- Orientation 2: N]chr:pos]  (nucleotide, close bracket, mate, close bracket) --

test_that("parse_bnd_alt handles N]chr:pos] orientation", {
  result <- ggmethylation:::parse_bnd_alt("T]chr5:500000]")
  expect_equal(result$mate_chrom, "chr5")
  expect_equal(result$mate_pos,   500000L)
})

test_that("parse_bnd_alt handles G]chr10:1] orientation", {
  result <- ggmethylation:::parse_bnd_alt("G]chr10:1]")
  expect_equal(result$mate_chrom, "chr10")
  expect_equal(result$mate_pos,   1L)
})

# --- Orientation 3: [chr:pos[N  (open bracket, mate, open bracket, nucleotide) --

test_that("parse_bnd_alt handles [chr:pos[N orientation", {
  result <- ggmethylation:::parse_bnd_alt("[chr7:200[C")
  expect_equal(result$mate_chrom, "chr7")
  expect_equal(result$mate_pos,   200L)
})

test_that("parse_bnd_alt handles [chrX:123456[T orientation", {
  result <- ggmethylation:::parse_bnd_alt("[chrX:123456[T")
  expect_equal(result$mate_chrom, "chrX")
  expect_equal(result$mate_pos,   123456L)
})

# --- Orientation 4: ]chr:pos]N  (close bracket, mate, close bracket, nucleotide) --

test_that("parse_bnd_alt handles ]chr:pos]N orientation", {
  result <- ggmethylation:::parse_bnd_alt("]chr1:1000000]G")
  expect_equal(result$mate_chrom, "chr1")
  expect_equal(result$mate_pos,   1000000L)
})

test_that("parse_bnd_alt handles ]chrY:500]A orientation", {
  result <- ggmethylation:::parse_bnd_alt("]chrY:500]A")
  expect_equal(result$mate_chrom, "chrY")
  expect_equal(result$mate_pos,   500L)
})

# --- Return types -------------------------------------------------------------

test_that("parse_bnd_alt mate_chrom is character", {
  result <- ggmethylation:::parse_bnd_alt("A[chr2:100[")
  expect_type(result$mate_chrom, "character")
})

test_that("parse_bnd_alt mate_pos is integer", {
  result <- ggmethylation:::parse_bnd_alt("A[chr2:100[")
  expect_type(result$mate_pos, "integer")
})

test_that("parse_bnd_alt returns a list with two named elements", {
  result <- ggmethylation:::parse_bnd_alt("A[chr2:100[")
  expect_true(is.list(result))
  expect_named(result, c("mate_chrom", "mate_pos"))
})

# --- Large chromosomal positions (no overflow) --------------------------------

test_that("parse_bnd_alt handles large positions correctly", {
  result <- ggmethylation:::parse_bnd_alt("N[chr1:248956422[")
  expect_equal(result$mate_chrom, "chr1")
  expect_equal(result$mate_pos,   248956422L)
})
