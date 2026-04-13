# Tests for build_variant_overlay() in R/variant_overlay.R
# Internal function — call via ggmethylation:::build_variant_overlay()
#
# Signature:
#   build_variant_overlay(data, variants, bnd_match_tol = 50L)
#
# Returns: named list(snv, sv, bnd, sa_reads), each element is a list of
#          ggplot2 layers or NULL.  Returns NULL when variants is NULL or
#          not a variant_data object.

# --- Test 1: NULL variants returns NULL ---

test_that("build_variant_overlay returns NULL when variants is NULL", {
  expect_null(ggmethylation:::build_variant_overlay(NULL, NULL))
})

# --- Test 2: Non-variant_data input returns NULL ---

test_that("build_variant_overlay returns NULL for non-variant_data input", {
  plain_list <- list(variants = data.frame(position = 1L))
  expect_null(ggmethylation:::build_variant_overlay(NULL, plain_list))
})

# --- Test 3: variant_data with 0-row variants returns list with all-NULL elements ---

test_that("build_variant_overlay returns named list with NULL elements for empty variant_data", {
  empty_vd <- structure(
    list(
      variants = data.frame(
        position    = integer(0),
        ref         = character(0),
        alt         = character(0),
        type        = character(0),
        end         = integer(0),
        mate_chrom  = character(0),
        mate_pos    = integer(0),
        stringsAsFactors = FALSE
      ),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
    ),
    class = "variant_data"
  )

  # Build a minimal methylation_data so the function doesn't error on data$reads etc.
  minimal_md <- structure(
    list(
      reads = data.frame(
        read_name = character(0), start = integer(0), end = integer(0),
        strand = character(0), lane = integer(0),
        stringsAsFactors = FALSE
      ),
      sites = data.frame(
        read_name = character(0), position = integer(0),
        mod_prob  = numeric(0), stringsAsFactors = FALSE
      ),
      sequences = NULL,
      cigars    = NULL,
      region    = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
    ),
    class = "methylation_data"
  )

  result <- ggmethylation:::build_variant_overlay(minimal_md, empty_vd)

  expect_equal(names(result), c("snv", "sv", "bnd", "sa_reads"))
  expect_null(result$snv)
  expect_null(result$sv)
  expect_null(result$bnd)
  expect_null(result$sa_reads)
})

# --- Test 4: variant_data with one SNV row + methylation_data without sequences
#             produces a warning and returns snv = NULL ---

test_that("build_variant_overlay warns and returns snv=NULL when sequences unavailable", {
  minimal_md <- structure(
    list(
      reads = data.frame(
        read_name = "r1", start = 100L, end = 200L,
        strand = "+", lane = 1L,
        stringsAsFactors = FALSE
      ),
      sites = data.frame(
        read_name = character(0), position = integer(0),
        mod_prob  = numeric(0), stringsAsFactors = FALSE
      ),
      sequences = NULL,   # forces sequences-unavailable path
      cigars    = NULL,
      region    = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 300))
    ),
    class = "methylation_data"
  )

  snv_vd <- structure(
    list(
      variants = data.frame(
        position   = 150L,
        ref        = "A",
        alt        = "T",
        type       = "SNV",
        end        = 150L,
        mate_chrom = NA_character_,
        mate_pos   = NA_integer_,
        stringsAsFactors = FALSE
      ),
      region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 300))
    ),
    class = "variant_data"
  )

  expect_warning(
    result <- ggmethylation:::build_variant_overlay(minimal_md, snv_vd)
  )
  expect_null(result$snv)
})
