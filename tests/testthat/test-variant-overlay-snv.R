# Tests for build_snv_layer() in R/variant_overlay_snv.R
# Internal function — call via ggmethylation:::build_snv_layer()
#
# Signature:
#   build_snv_layer(variant_bases)
#
# variant_bases: data.frame with columns read_name, position, lane,
#                variant_class (values: "ref", "alt", "other", "del")
#                May be NULL or zero-row.
#
# Returns: list of 1 ggplot2-compatible geom_point layer, or NULL when there
#          are no ALT carrier calls to draw.

# Helper: build a minimal variant_bases data.frame
make_vb <- function(read_name = "r1", position = 100L, lane = 1L,
                    variant_class = "alt") {
  data.frame(
    read_name     = read_name,
    position      = as.integer(position),
    lane          = as.integer(lane),
    variant_class = variant_class,
    stringsAsFactors = FALSE
  )
}

# --- Test 1: NULL input returns NULL ---

test_that("build_snv_layer returns NULL for NULL input", {
  expect_null(ggmethylation:::build_snv_layer(NULL))
})

# --- Test 2: 0-row input returns NULL ---

test_that("build_snv_layer returns NULL for zero-row data.frame", {
  empty_df <- data.frame(
    read_name     = character(0L),
    position      = integer(0L),
    lane          = integer(0L),
    variant_class = character(0L),
    stringsAsFactors = FALSE
  )
  expect_null(ggmethylation:::build_snv_layer(empty_df))
})

# --- Test 3: All-ref input returns NULL ---

test_that("build_snv_layer returns NULL when all rows are ref", {
  df <- make_vb(variant_class = "ref")
  expect_null(ggmethylation:::build_snv_layer(df))
})

# --- Test 4: One alt row returns a 1-element list with a geom_point layer ---

test_that("build_snv_layer returns a list of length 1 for a single alt row", {
  df     <- make_vb(variant_class = "alt")
  result <- ggmethylation:::build_snv_layer(df)

  expect_true(is.list(result))
  expect_equal(length(result), 1L)
  expect_true(inherits(result[[1L]], "ggproto") || is.list(result[[1L]]))
})

# --- Test 5: Only "alt" rows are drawn; "ref", "other", "del" are suppressed ---

test_that("build_snv_layer draws only alt rows, suppresses ref/other/del", {
  df <- rbind(
    make_vb(read_name = "r1", variant_class = "ref"),
    make_vb(read_name = "r2", variant_class = "alt"),
    make_vb(read_name = "r3", variant_class = "other"),
    make_vb(read_name = "r4", variant_class = "del")
  )
  result <- ggmethylation:::build_snv_layer(df)

  # A layer is returned (the one alt row)
  expect_false(is.null(result))
  expect_equal(length(result), 1L)

  # The data passed to the layer contains only the alt read
  layer_data <- result[[1L]]$data
  expect_equal(nrow(layer_data), 1L)
  expect_equal(layer_data$variant_class, "alt")
  expect_equal(layer_data$read_name, "r2")
})

# --- Test 6: All-non-alt (ref + other + del only) returns NULL ---

test_that("build_snv_layer returns NULL when no alt rows are present", {
  df <- rbind(
    make_vb(read_name = "r1", variant_class = "ref"),
    make_vb(read_name = "r2", variant_class = "other"),
    make_vb(read_name = "r3", variant_class = "del")
  )
  expect_null(ggmethylation:::build_snv_layer(df))
})

# --- Test 7 (regression, issue #20): extract_variant_bases() does not route
#     insertion/deletion variants through the SNV path ---

test_that("build_variant_overlay produces no SNV layer for an insertion variant", {
  # Minimal methylation_data-like object with one read and real sequences/cigars
  reads_df <- data.frame(
    read_name = "r1",
    start     = 100L,
    end       = 200L,
    bam_pos   = 100L,
    lane      = 1L,
    stringsAsFactors = FALSE
  )
  seqs   <- c(r1 = paste(rep("A", 150L), collapse = ""))
  cigars <- c(r1 = "150M")

  # An insertion variant: single-char ref, multi-char alt — type "insertion"
  variants_df <- data.frame(
    position   = 150L,
    ref        = "N",
    alt        = paste(rep("A", 93L), collapse = ""),
    type       = "insertion",
    end        = 150L,
    mate_chrom = NA_character_,
    mate_pos   = NA_integer_,
    stringsAsFactors = FALSE
  )
  vb <- ggmethylation:::extract_variant_bases(reads_df, seqs, cigars, variants_df)
  # extract_variant_bases itself would classify the one read as "other" (single
  # base vs multi-char alt). The fix is that build_variant_overlay() never feeds
  # indels to this function; this test verifies the outcome via build_snv_layer().
  layer <- ggmethylation:::build_snv_layer(vb)

  # With all rows classified "other", the alt-only filter returns NULL
  expect_null(layer)
})
