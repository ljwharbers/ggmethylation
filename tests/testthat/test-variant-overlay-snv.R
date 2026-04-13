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
# Returns: list of 4 ggplot2-compatible objects, or NULL when there are no
#          non-ref calls to draw.

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

# --- Test 4: One alt row returns a 4-element list of ggplot2-compatible objects ---

test_that("build_snv_layer returns a list of length 4 for a single alt row", {
  df     <- make_vb(variant_class = "alt")
  result <- ggmethylation:::build_snv_layer(df)

  expect_true(is.list(result))
  expect_equal(length(result), 4L)

  # Each element should be a ggplot2-compatible object:
  #   - "ggproto"  for layer / scale objects
  #   - "new_aes"  for ggnewscale::new_scale_colour() (a classed character)
  #   - list       fallback for future ggplot2 internals
  for (elem in result) {
    expect_true(
      inherits(elem, "ggproto") || is.list(elem) || inherits(elem, "new_aes"),
      info = paste("Element class:", paste(class(elem), collapse = ", "))
    )
  }
})

# --- Test 5: Mixed input (ref + alt + other) — ref rows filtered out, list returned ---

test_that("build_snv_layer filters ref rows and returns a list for mixed input", {
  df <- rbind(
    make_vb(read_name = "r1", variant_class = "ref"),
    make_vb(read_name = "r2", variant_class = "alt"),
    make_vb(read_name = "r3", variant_class = "other")
  )
  result <- ggmethylation:::build_snv_layer(df)

  expect_false(is.null(result))
  expect_true(is.list(result))
  expect_equal(length(result), 4L)
})
