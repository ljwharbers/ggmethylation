# Tests for match_sa_to_vcf_bnd() in R/variant_overlay_bnd.R
# Internal function — call via ggmethylation:::match_sa_to_vcf_bnd()
#
# Signature:
#   match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)
#
# reads:  data.frame with columns read_name, sa_chrom, sa_pos (and optionally
#         others). sa_chrom / sa_pos are NA for non-SA reads.
# bnd_df: data.frame with columns position (int), mate_chrom (chr),
#         mate_pos (int). May be NULL or zero-row.
# tol:    integer position tolerance in bp (default 50L).
#
# Returns: reads with new logical column vcf_validated.

# Helpers --------------------------------------------------------------------

make_reads <- function(read_name = "r1",
                       sa_chrom  = "chr7",
                       sa_pos    = 1000L) {
  data.frame(
    read_name = read_name,
    sa_chrom  = sa_chrom,
    sa_pos    = as.integer(sa_pos),
    stringsAsFactors = FALSE
  )
}

make_bnd <- function(position   = 1010L,
                     mate_chrom = "chr7",
                     mate_pos   = 999L) {
  data.frame(
    position   = as.integer(position),
    mate_chrom = mate_chrom,
    mate_pos   = as.integer(mate_pos),
    stringsAsFactors = FALSE
  )
}

# --- Test 1: NULL bnd_df — all reads get vcf_validated = FALSE ---------------

test_that("match_sa_to_vcf_bnd returns vcf_validated FALSE for NULL bnd_df", {
  reads  <- make_reads()
  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, NULL)

  expect_true("vcf_validated" %in% names(result))
  expect_false(any(result$vcf_validated))
})

# --- Test 2: 0-row bnd_df — all reads get vcf_validated = FALSE --------------

test_that("match_sa_to_vcf_bnd returns vcf_validated FALSE for 0-row bnd_df", {
  reads  <- make_reads()
  bnd_df <- data.frame(
    position   = integer(0L),
    mate_chrom = character(0L),
    mate_pos   = integer(0L),
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df)

  expect_true("vcf_validated" %in% names(result))
  expect_false(any(result$vcf_validated))
})

# --- Test 3: Matching — sa_chrom == mate_chrom AND |sa_pos - position| <= tol → TRUE ---

test_that("match_sa_to_vcf_bnd sets vcf_validated TRUE for a full match", {
  reads  <- make_reads(sa_chrom = "chr7",  sa_pos  = 1000L)
  bnd_df <- make_bnd( mate_chrom = "chr7", position = 1010L)  # delta = 10 <= 50

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_true(result$vcf_validated[1L])
})

# --- Test 4: Position-only match rejected — sa_chrom != mate_chrom → FALSE --

test_that("match_sa_to_vcf_bnd returns FALSE when chrom does not match", {
  reads  <- make_reads(sa_chrom = "chr8",  sa_pos  = 1000L)
  bnd_df <- make_bnd( mate_chrom = "chr7", position = 1010L)  # delta = 10, but wrong chrom

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

# --- Test 5: Chromosome-only match rejected — |sa_pos - position| > tol → FALSE ---

test_that("match_sa_to_vcf_bnd returns FALSE when position is out of tolerance", {
  reads  <- make_reads(sa_chrom = "chr7",  sa_pos  = 1000L)
  bnd_df <- make_bnd( mate_chrom = "chr7", position = 1100L)  # delta = 100 > 50

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

# --- Test 6: Non-SA reads (is.na(sa_chrom)) are always FALSE -----------------

test_that("match_sa_to_vcf_bnd returns FALSE for reads with no SA tag", {
  reads <- data.frame(
    read_name = c("r1", "r2"),
    sa_chrom  = c(NA_character_, NA_character_),
    sa_pos    = c(NA_integer_,   NA_integer_),
    stringsAsFactors = FALSE
  )
  bnd_df <- make_bnd(mate_chrom = "chr7", position = 1010L)

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(any(result$vcf_validated))
})

# --- Test 7: Exact boundary conditions — tol and tol+1 ----------------------

test_that("match_sa_to_vcf_bnd handles exact boundary: |delta| == tol → TRUE", {
  reads  <- make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df <- make_bnd( mate_chrom = "chr7", position = 1050L)  # delta = 50 == tol

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_true(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd handles exact boundary: |delta| == tol + 1 → FALSE", {
  reads  <- make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df <- make_bnd( mate_chrom = "chr7", position = 1051L)  # delta = 51 > tol

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

# --- Test 8: vcf_validated column is logical type ----------------------------

test_that("match_sa_to_vcf_bnd always produces a logical vcf_validated column", {
  reads  <- make_reads()
  bnd_df <- make_bnd()

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df)

  expect_type(result$vcf_validated, "logical")
})

# --- Tests for build_bnd_layer() ---------------------------------------------

test_that("build_bnd_layer returns NULL for NULL input", {
  expect_null(ggmethylation:::build_bnd_layer(NULL))
})

test_that("build_bnd_layer returns a list of length 2 for a one-row bnd_df", {
  bnd_df <- data.frame(
    position   = 1000L,
    mate_chrom = "chr7",
    mate_pos   = 5000L,
    stringsAsFactors = FALSE
  )
  result <- ggmethylation:::build_bnd_layer(bnd_df)

  expect_type(result, "list")
  expect_length(result, 2L)
})

# --- Bonus: mixed reads (some SA, some not) -----------------------------------

test_that("match_sa_to_vcf_bnd handles mixed SA and non-SA reads correctly", {
  reads <- data.frame(
    read_name = c("r1", "r2", "r3"),
    sa_chrom  = c("chr7", NA_character_, "chr7"),
    sa_pos    = c(1000L,   NA_integer_,   2000L),
    stringsAsFactors = FALSE
  )
  # Only r1 should match (correct chrom + within tol); r2 has no SA; r3 is out of tol
  bnd_df <- make_bnd(mate_chrom = "chr7", position = 1010L)

  result <- ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_true(result$vcf_validated[1L])   # r1: match
  expect_false(result$vcf_validated[2L])  # r2: no SA
  expect_false(result$vcf_validated[3L])  # r3: wrong position
})
