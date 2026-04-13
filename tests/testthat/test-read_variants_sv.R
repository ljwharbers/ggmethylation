# Tests for symbolic SV classification in R/read_variants.R
# These test classify_variant_row() directly — no VCF file required.

# --- <DEL> -------------------------------------------------------------------

test_that("classify_variant_row gives type DEL for <DEL> ALT", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L
  )
  expect_equal(result$type, "DEL")
})

test_that("classify_variant_row uses INFO/END for DEL end coordinate", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L, end = 5000L
  )
  expect_equal(result$end, 5000L)
})

test_that("classify_variant_row falls back to pos + abs(SVLEN) when END absent", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L, end = NA_integer_, svlen = -4000
  )
  expect_equal(result$end, 5000L)  # 1000 + abs(-4000)
})

test_that("classify_variant_row uses abs(SVLEN) so negative SVLEN is handled", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 2000L, end = NA_integer_, svlen = -1000
  )
  expect_equal(result$end, 3000L)  # 2000 + abs(-1000)
})

test_that("classify_variant_row falls back to pos when both END and SVLEN absent", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L
  )
  expect_equal(result$end, 1000L)
})

test_that("classify_variant_row DEL has NA mate columns", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L, end = 5000L
  )
  expect_equal(result$mate_chrom, NA_character_)
  expect_equal(result$mate_pos,   NA_integer_)
})

# --- INFO/END takes priority over SVLEN --------------------------------------

test_that("classify_variant_row prefers END over SVLEN when both present", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L, end = 2000L, svlen = -9999
  )
  expect_equal(result$end, 2000L)
})

# --- <DUP> -------------------------------------------------------------------

test_that("classify_variant_row gives type DUP for <DUP> ALT", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DUP>", ref = "N", pos = 500L, end = 1500L
  )
  expect_equal(result$type, "DUP")
  expect_equal(result$end,  1500L)
})

test_that("classify_variant_row DUP SVLEN fallback works", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DUP>", ref = "N", pos = 500L, svlen = 1000
  )
  expect_equal(result$type, "DUP")
  expect_equal(result$end,  1500L)
})

# --- <INV> -------------------------------------------------------------------

test_that("classify_variant_row gives type INV for <INV> ALT", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<INV>", ref = "N", pos = 300L, end = 800L
  )
  expect_equal(result$type, "INV")
  expect_equal(result$end,  800L)
})

# --- SVLEN supplied as a list (VariantAnnotation list-type INFO field) -------

test_that("classify_variant_row handles list-type SVLEN correctly", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L, end = NA_integer_,
    svlen = list(-3000)
  )
  expect_equal(result$end, 4000L)  # 1000 + abs(-3000)
})

# --- Return types for SV records ---------------------------------------------

test_that("classify_variant_row SV end is integer", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DEL>", ref = "N", pos = 1000L, end = 5000L
  )
  expect_type(result$end, "integer")
})

test_that("classify_variant_row SV type is character", {
  result <- ggmethylation:::classify_variant_row(
    alt = "<DUP>", ref = "N", pos = 100L, end = 200L
  )
  expect_type(result$type, "character")
})

# --- Standard types are unaffected -------------------------------------------

test_that("classify_variant_row SNV classification preserved", {
  result <- ggmethylation:::classify_variant_row(alt = "A", ref = "T", pos = 100L)
  expect_equal(result$type, "SNV")
  expect_equal(result$end,  100L)
})

test_that("classify_variant_row insertion classification preserved", {
  result <- ggmethylation:::classify_variant_row(alt = "ATCG", ref = "A", pos = 200L)
  expect_equal(result$type, "insertion")
  expect_equal(result$end,  200L)
})

test_that("classify_variant_row deletion classification preserved", {
  result <- ggmethylation:::classify_variant_row(alt = "A", ref = "ATCG", pos = 300L)
  expect_equal(result$type, "deletion")
  expect_equal(result$end,  300L)
})

# --- BND via SVTYPE field ----------------------------------------------------

test_that("classify_variant_row classifies BND when svtype is BND", {
  result <- ggmethylation:::classify_variant_row(
    alt = "N", ref = "N", pos = 500L, svtype = "BND"
  )
  expect_equal(result$type, "BND")
  expect_equal(result$end,  500L)
})

test_that("classify_variant_row classifies BND from bracket ALT even without svtype", {
  result <- ggmethylation:::classify_variant_row(
    alt = "A[chr2:12345[", ref = "A", pos = 100L
  )
  expect_equal(result$type,       "BND")
  expect_equal(result$mate_chrom, "chr2")
  expect_equal(result$mate_pos,   12345L)
})

# --- NA_character_ alt/ref guard ---------------------------------------------

test_that("classify_variant_row returns NA_character_ type for NA alt", {
  result <- ggmethylation:::classify_variant_row(
    alt    = NA_character_,
    ref    = "A",
    pos    = 100L,
    svtype = NA_character_,
    end    = NA,
    svlen  = NA
  )
  expect_equal(result$type,       NA_character_)
  expect_equal(result$end,        100L)
  expect_equal(result$mate_chrom, NA_character_)
  expect_equal(result$mate_pos,   NA_integer_)
})
