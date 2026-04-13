# Tests for build_sv_layer() in R/variant_overlay_sv.R
# Internal function — call via ggmethylation:::build_sv_layer()
#
# Signature:
#   build_sv_layer(reads, sv_df, region_start, region_end)
#
# reads:        data.frame with columns read_name, start, end, lane
# sv_df:        data.frame with columns position (int), end (int), type (chr)
#               May be NULL or zero-row.
# region_start: integer — left boundary of plot window
# region_end:   integer — right boundary of plot window
#
# Returns: flat list of ggplot2 layer objects, or NULL when no overlaps.

# Helpers --------------------------------------------------------------------

make_reads <- function(read_name = "r1", start = 100L, end = 400L, lane = 1L) {
  data.frame(
    read_name = read_name,
    start     = as.integer(start),
    end       = as.integer(end),
    lane      = as.integer(lane),
    stringsAsFactors = FALSE
  )
}

make_sv <- function(position = 200L, end = 300L, type = "DEL") {
  data.frame(
    position = as.integer(position),
    end      = as.integer(end),
    type     = type,
    stringsAsFactors = FALSE
  )
}

# --- Test 1: NULL sv_df returns NULL ----------------------------------------

test_that("build_sv_layer returns NULL for NULL sv_df", {
  reads <- make_reads()
  expect_null(ggmethylation:::build_sv_layer(reads, NULL, 100L, 500L))
})

# --- Test 2: 0-row sv_df returns NULL ----------------------------------------

test_that("build_sv_layer returns NULL for zero-row sv_df", {
  reads  <- make_reads()
  sv_df <- data.frame(
    position = integer(0L),
    end      = integer(0L),
    type     = character(0L),
    stringsAsFactors = FALSE
  )
  expect_null(ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L))
})

# --- Test 3: No overlapping reads returns NULL ------------------------------

test_that("build_sv_layer returns NULL when reads do not overlap the SV", {
  # Read spans 100-200; SV spans 300-400 — no overlap
  reads <- make_reads(start = 100L, end = 200L)
  sv_df <- make_sv(position = 300L, end = 400L)
  expect_null(ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L))
})

# --- Test 4: Basic DEL -------------------------------------------------------

test_that("build_sv_layer returns a non-NULL list for a basic DEL", {
  reads  <- make_reads(start = 100L, end = 400L)
  sv_df  <- make_sv(position = 200L, end = 300L, type = "DEL")
  result <- ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L)

  expect_false(is.null(result))
  expect_true(is.list(result))
  expect_gte(length(result), 1L)
})

# --- Test 5: Clip at left boundary -------------------------------------------

test_that("build_sv_layer handles clipped_left when SV position < region_start", {
  # SV starts at 50 (< region_start 100), ends at 300
  reads  <- make_reads(start = 100L, end = 400L)
  sv_df  <- make_sv(position = 50L, end = 300L, type = "DEL")
  result <- ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L)

  expect_false(is.null(result))
  expect_true(is.list(result))
})

# --- Test 6: Clip at right boundary ------------------------------------------

test_that("build_sv_layer handles clipped_right when SV end > region_end", {
  # SV ends at 600 (> region_end 500), starts at 200
  reads  <- make_reads(start = 200L, end = 600L)
  sv_df  <- make_sv(position = 200L, end = 600L, type = "DEL")
  result <- ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L)

  expect_false(is.null(result))
  expect_true(is.list(result))
})

# --- Test 7: DUP type returns non-NULL list ----------------------------------

test_that("build_sv_layer returns layers for DUP type", {
  reads  <- make_reads(start = 100L, end = 500L)
  sv_df  <- make_sv(position = 200L, end = 400L, type = "DUP")
  result <- ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L)

  expect_false(is.null(result))
  expect_true(is.list(result))
})

# --- Test 8: INV type returns non-NULL list ----------------------------------

test_that("build_sv_layer returns layers for INV type", {
  reads  <- make_reads(start = 100L, end = 500L)
  sv_df  <- make_sv(position = 150L, end = 350L, type = "INV")
  result <- ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L)

  expect_false(is.null(result))
  expect_true(is.list(result))
})

# --- Test 9: Multiple reads, multiple SVs -----------------------------------

test_that("build_sv_layer handles multiple reads and multiple SVs", {
  reads <- make_reads(
    read_name = c("r1", "r2"),
    start     = c(100L, 300L),
    end       = c(400L, 700L),
    lane      = c(1L, 2L)
  )
  sv_df <- rbind(
    make_sv(position = 150L, end = 250L, type = "DEL"),
    make_sv(position = 350L, end = 600L, type = "DUP")
  )
  result <- ggmethylation:::build_sv_layer(reads, sv_df, 100L, 800L)

  expect_false(is.null(result))
  expect_true(is.list(result))
})

# --- Test 10: SV entirely outside plot window returns NULL ------------------

test_that("build_sv_layer returns NULL when SV is fully outside the plot window", {
  reads <- make_reads(start = 100L, end = 400L)
  # SV at 600-700 is outside region 100-500
  sv_df <- make_sv(position = 600L, end = 700L, type = "INV")
  expect_null(ggmethylation:::build_sv_layer(reads, sv_df, 100L, 500L))
})
